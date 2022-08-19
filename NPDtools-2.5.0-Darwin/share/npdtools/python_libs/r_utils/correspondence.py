import time
import config
from os.path import join, dirname, basename, isdir
from log_utils import info, error, warning
from common import check_file
from fname_utils import _get_corresp_fpath, compare_basename_wo_extension
from collections import OrderedDict
from ordered_set import OrderedSet

# for downloading from NCBI
import socket
socket.setdefaulttimeout(10.0)
try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen
import xml.etree.ElementTree as ET
from common import warning

connection_errors = 0


def get_spec_seq_correspondence(cfg, corr_fpath, spectra_fpaths, sequence_fpaths, file_mapping=None):
    spec_keyword = 'spectr'
    seq_keyword = 'sequence'
    spec_idx = 0
    seq_idx = 1

    def __get_spec_sec_idx(header_line):
        global spec_idx
        global seq_idx
        if seq_keyword in header_line.lower() and spec_keyword in header_line.lower():  # header is present
            header_columns = header_line.strip().lower().split('\t')
            for idx, column in enumerate(header_columns):
                if spec_keyword in column.lower():
                    spec_idx = idx
                elif seq_keyword in column.lower():
                    seq_idx = idx
            return True
        return False

    def __get_real_fpath(fname, real_fpaths, spectra=True):
        if spectra:
            exts = config.allowed_spectrum_extensions
        else:
            exts = config.allowed_sequence_extensions

        for real_fpath in real_fpaths:
            fpath_to_compare = basename(real_fpath)
            if file_mapping and fpath_to_compare in file_mapping:
                fpath_to_compare = basename(file_mapping[fpath_to_compare])
            if compare_basename_wo_extension(fpath_to_compare, fname, exts):
                return real_fpath
        return None

    def __get_dict_of_sets(spectra_fpaths, sequence_fpaths):
        dict_of_sets = OrderedDict()
        for spectra_fpath in spectra_fpaths:
            dict_of_sets[spectra_fpath] = OrderedSet(sequence_fpaths)
        return dict_of_sets

    if corr_fpath and check_file(corr_fpath):
        spec_seq_corr = __get_dict_of_sets(spectra_fpaths, [])

        with open(corr_fpath) as f:
            content = f.readlines()
        if __get_spec_sec_idx(content[0]):
            del content[0]  # removing header

        for spec_seq_corr_line in content:
            if not spec_seq_corr_line.strip():
                continue
            spec_seq_corr_line = spec_seq_corr_line.strip()
            if len(spec_seq_corr_line.split('\t')) <= max(spec_idx, seq_idx):
                warning('Incorrect entry in correspondence file ({corr_fpath})! '
                        'Probably spaces were used instead of tabs! '
                        'Skipping: {spec_seq_corr_line}'.format(**locals()))
                continue

            spec_fname = spec_seq_corr_line.split('\t')[spec_idx].strip()
            seq_fname = spec_seq_corr_line.split('\t')[seq_idx].strip()

            spectra_fpath = __get_real_fpath(spec_fname, spectra_fpaths)
            if spectra_fpath is None:
                warning("Spectra filename ({spec_fname}) from correspondence file ({corr_fpath}) "
                        "was not found among input spectra files! Skipping this entry: {spec_seq_corr_line}".format(**locals()))
                continue
            if seq_fname.startswith('#') and ':' in seq_fname:  # special case
                info('Trying to download sequence file ({seq_fname}) specified in correspondence file ({corr_fpath}) '
                     'from NCBI'.format(**locals()))
                sequence_fpath = join(cfg.work_dirpath, seq_fname.split(':')[1] + '.fasta')
                if check_file(sequence_fpath):
                    info('  Sequence file was already downloaded, reusing ({sequence_fpath})!'.format(**locals()))
                else:
                    if download_by_id(seq_fname, sequence_fpath) is None:
                        warning('  Failed to download sequence file! '
                                'Skipping this entry: {spec_seq_corr_line}'.format(**locals()))
                        continue
                    else:
                        info('  Successfully downloaded and saved to {sequence_fpath}'.format(**locals()))
            else:
                sequence_fpath = __get_real_fpath(seq_fname, sequence_fpaths, spectra=False)
                if sequence_fpath is None:
                    warning('Sequence filename ({seq_fname}) from correspondence file ({corr_fpath}) '
                            'was not found among input sequence files! Skipping this entry: {spec_seq_corr_line}'.format(**locals()))
                    continue
            spec_seq_corr[spectra_fpath].add(sequence_fpath)
    else:
        spec_seq_corr = __get_dict_of_sets(spectra_fpaths, sequence_fpaths)

    return spec_seq_corr


def save_spec_seq_correspondence(cfg, spec_seq_corr):
    corresp_fpath = _get_corresp_fpath(cfg)
    with open(corresp_fpath, 'w') as f:
        f.write('SpectraFile\tSequenceFile\n')
        for spec_fpath, set_of_seq in spec_seq_corr.items():
            for seq_fpath in set_of_seq:
                f.write(spec_fpath + '\t' + seq_fpath + '\n')
    return corresp_fpath


def try_send_request(url):
    attempts = 0
    response = None
    global connection_errors
    while attempts < 3:
        try:
            request = urlopen(url)
            connection_errors = 0
            response = request.read()
            if not isinstance(response, str):
                response = response.decode('utf-8')
            if response is None or 'ERROR' in response:
                request.close()
                raise Exception
            break
        except Exception:
            # _, exc_value, _ = sys.exc_info()
            # logger.exception(exc_value)
            attempts += 1
            if attempts >= 3:
                connection_errors += 1
                if connection_errors >= 3:
                    warning('Cannot established internet connection to download reference genomes!')
                return None
            # NCBI recommends users post no more than three URL requests per second, so adding artificial 1-sec delay
            # see more in similar Quast issue: https://github.com/ablab/quast/issues/8
            time.sleep(1)
    return response


def download_by_id(ref_id, output_fpath):
    ncbi_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    type_, id_ = ref_id.split(':')
    response = try_send_request(ncbi_url + 'esearch.fcgi?db=assembly&term=' + id_)
    xml_tree = ET.fromstring(response)
    if xml_tree.find('Count').text == '0':  # Organism is not found
        return None
    assembly_id = xml_tree.find('IdList').find('Id')
    if assembly_id is None:
        return None
    if 'GenBank' in type_:
        link_name = 'assembly_nuccore_insdc'
    else:
        link_name = 'assembly_nuccore_refseq'
    assembly_id = assembly_id.text
    response = try_send_request(
        ncbi_url + 'elink.fcgi?dbfrom=assembly&db=nuccore&id=' + assembly_id + '&linkname="' + link_name + '"')
    xml_tree = ET.fromstring(response)

    link_set = xml_tree.find('LinkSet')
    if link_set is None:
        return None

    link_db = xml_tree.find('LinkSet').find('LinkSetDb')
    if link_db is None:
        return None

    fasta_ids = [link.find('Id').text for link in link_db.findall('Link')]
    fasta_files = []
    for fasta_id in fasta_ids:
        fasta = try_send_request(ncbi_url + 'efetch.fcgi?db=sequences&id=%s&rettype=fasta&retmode=text' % fasta_id)
        if fasta:
            fasta_files.append(fasta)
    is_first_piece = False
    with open(output_fpath, "w") as fasta_file:
        for fasta in fasta_files:
            if not is_first_piece:
                is_first_piece = True
            else:
                fasta = '\n' + fasta.rstrip()
            fasta_file.write(fasta.rstrip())
    return output_fpath