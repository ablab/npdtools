#!/usr/bin/env python

import collections
from glob import glob
from os.path import join, dirname, splitext, basename, isdir, isfile

from log_utils import info, error, warning
from fname_utils import _get_fpath_by_idx
from file_utils import verify_file, verify_dir
from common import parse_params_xml
from r_utils.pipeline import PSM_ORF
import config


def check_spec_network_files(basedir):
    verify_dir(basedir, 'SpecNets root')
    params, nodes, edges = get_spec_network_files(basedir)
    verify_file(params, 'SpecNets params')
    verify_file(nodes, 'SpecNets cluster info')
    verify_file(edges, 'SpecNets pairs info')


def get_spec_network_files(basedir):
    params = join(basedir, 'params.xml')
    clusters_dir = join(basedir, 'clusterinfo')
    verify_dir(clusters_dir, 'SpecNets clusters')
    clusters = glob(join(clusters_dir, '*.clusterinfo'))
    if len(clusters) != 1:
        error('SpecNets clusters dir (%s) should contain exactly one *.clusterinfo file!' % clusters_dir)
    edges_dir = join(basedir, 'networkedges_selfloop')
    verify_dir(edges_dir, 'SpecNets edges')
    edges = glob(join(edges_dir, '*.pairsinfo'))
    if not len(edges):
        edges = glob(join(edges_dir, '*.selfloop'))  # format for new METABOLOMICS-SNETS-V2 workflow
    if len(edges) != 1:
        error('SpecNets edges dir (%s) should contain exactly one *.pairsinfo or *.selfloop file!' % edges_dir)
    return params, clusters[0], edges[0]


def __get_spec_id(psm):
    spec_key = basename(psm.spectrum_fpath)
    scan_key = psm.scan_id
    return spec_key, scan_key


# utility for Dereplicator
def parse_and_propagate(cfg, spectra_fpaths, spec_nets_dir, all_PSMs):
    info("\n\tStage 2: parsing SpecNets and stage 1 results to create list of putative PSMs")
    params, nodes, edges = get_spec_network_files(spec_nets_dir)
    # parsing SpecNets
    _, spec_nets_file_mapping = parse_params_xml(params)
    cluster_to_specs = get_cluster_spec_corresp(nodes, spec_nets_file_mapping)
    edges = get_edges(edges)
    # parsing Dereplicator
    dereplicator_spec_to_compounds = get_spec_compound_corresp(
        [psm for psm in all_PSMs if float(psm.p_value) < config.dereplicator_p_value_threshold])
    # propagate
    spec_putative_compounds = propagate_annotations(dereplicator_spec_to_compounds, cluster_to_specs, edges)
    write_putative_compounds(cfg, spectra_fpaths, spec_putative_compounds)


# utility for MetaMiner
def parse_and_report_propagation(cfg, spectra_fpaths, spec_nets_dir, important_PSMs):
    info("\n\tPostprocessing: parsing SpecNets and identification results to propagate them via the network")
    params, nodes, edges = get_spec_network_files(spec_nets_dir)
    # parsing SpecNets
    _, spec_nets_file_mapping = parse_params_xml(params)
    cluster_to_spec, spec_to_metadata = get_cluster_spec_corresp(nodes, spec_nets_file_mapping, with_metadata=True)
    edges, edge_to_metadata = get_edges(edges, with_metadata=True)
    spec_to_cluster = get_spec_cluster_corresp(cluster_to_spec)
    cluster_to_neighbours = get_clusters_adjacency_list(edges)
    # preprocessing
    # filtering clusters (joining isotopic shifts and primary spectra)
    merge_isotopes(cluster_to_spec, spec_to_cluster, edges, cluster_to_neighbours)
    # filtering PSMs (only the best PSM per spectrum)
    important_PSMs = get_unique_PSMs_by_spectrum(important_PSMs)

    # propagate
    results_fpaths = write_propagations(join(cfg.output_dirpath, 'propagations'), important_PSMs,
                                        spec_to_cluster, cluster_to_spec, cluster_to_neighbours, edges,
                                        spec_to_metadata, edge_to_metadata)
    results_fpaths_string = ", ".join([results_fpaths[0]] + list(map(basename, results_fpaths[1:])))
    info("\n\tPostprocessing done: results are in " + results_fpaths_string)


def merge_isotopes(cluster_to_spec, spec_to_cluster, edges, cluster_to_neighbours, max_delta=2.5, min_cosine=0.0):
    # GENERAL RULE: join two clusters if their abs(delta) <= max_delta and cosine >= min_cosine.
    # The new cluster is the lightest one out of two
    #
    # spec_to_cluster[spec_key][scan] = (cluster, info)             # info = (pm, rt)
    #
    # cluster_to_neighbours[cluster1Id].append((cluster2Id, info))  # info = (delta, cosine)
    # cluster_to_neighbours[cluster2Id].append((cluster1Id, info))
    # TODO: implement something
    pass


def get_unique_PSMs_by_spectrum(PSMs):
    PSMs.sort(key=lambda x: x.p_value)
    unique_PSMs = []
    unique_spec_scan_pairs = set()
    for psm in PSMs:
        spec_id = __get_spec_id(psm)
        if spec_id not in unique_spec_scan_pairs:
            unique_spec_scan_pairs.add(spec_id)
            unique_PSMs.append(psm)
    return unique_PSMs


def get_spec_cluster_corresp(cluster_to_spec):
    spec_to_cluster = dict()
    for cluster, specs in cluster_to_spec.items():
        for spec_id in specs:
            spec_to_cluster[spec_id] = cluster
    return spec_to_cluster


def get_clusters_adjacency_list(edges):
    cluster_to_neighbours = dict()
    for cluster1, cluster2 in edges:
        if cluster1 not in cluster_to_neighbours:
            cluster_to_neighbours[cluster1] = []
        if cluster2 not in cluster_to_neighbours:
            cluster_to_neighbours[cluster2] = []
        cluster_to_neighbours[cluster1].append(cluster2)
        cluster_to_neighbours[cluster2].append(cluster1)
    return cluster_to_neighbours


def get_spec_compound_corresp(psms):
    spec_id_to_compound_id = dict()
    for psm in psms:
        spec_id = __get_spec_id(psm)
        value = psm.compound_idx
        if spec_id not in spec_id_to_compound_id:
            spec_id_to_compound_id[spec_id] = set()
        spec_id_to_compound_id[spec_id].add(value)
    return spec_id_to_compound_id


def get_cluster_spec_corresp(nodes_fpath, spec_nets_file_mapping, with_metadata=False):
    cluster_id_to_spec_ids = dict()
    if with_metadata:
        spec_to_metadata = dict()
    with open(nodes_fpath) as f:
        header_clns = f.readline().split()
        clusterId_cln = header_clns.index('#ClusterIdx')
        fname_cln = header_clns.index('#Filename')
        scan_cln = header_clns.index('#Scan')
        pm_cln = header_clns.index('#ParentMass')
        rt_cln = header_clns.index('#RetTime')
        for line in f:
            clns = line.split()
            cluster_id = clns[clusterId_cln]
            spec_id = (basename(spec_nets_file_mapping[basename(clns[fname_cln])]), clns[scan_cln])
            if cluster_id not in cluster_id_to_spec_ids:
                cluster_id_to_spec_ids[cluster_id] = list()
            if with_metadata:
                pm = float(clns[pm_cln])
                rt = float(clns[rt_cln])
                spec_to_metadata[spec_id] = (pm, rt)
            cluster_id_to_spec_ids[cluster_id].append(spec_id)
    if with_metadata:
        return cluster_id_to_spec_ids, spec_to_metadata
    return cluster_id_to_spec_ids


def get_edges(edges_fpath, with_metadata=False):
    edges = []
    if with_metadata:
        edge_to_metadata = dict()
    with open(edges_fpath) as f:
        header_clns = f.readline().split()
        cluster1_cln = header_clns.index('CLUSTERID1')
        cluster2_cln = header_clns.index('CLUSTERID2')
        delta_cln = header_clns.index('DeltaMZ')
        cosine_cln = header_clns.index('Cosine')
        for line in f:
            clns = line.split()
            cluster1Id = clns[cluster1_cln]
            cluster2Id = clns[cluster2_cln]
            if cluster1Id == cluster2Id:  # selfloops section
                break
            if with_metadata:
                delta = float(clns[delta_cln])
                cosine = float(clns[cosine_cln])
                edge_to_metadata[min(cluster1Id, cluster2Id), max(cluster1Id, cluster2Id)] = (delta, cosine)
            edges.append((cluster1Id, cluster2Id))
    if with_metadata:
        return edges, edge_to_metadata
    return edges


def propagate_annotations(dereplicator_spec_to_compounds, cluster_to_specs, edges):
    cluster_to_compounds = dict()
    spec_putative_compounds = dict()
    for cluster, specs in cluster_to_specs.items():
        clust_compounds = set()
        for spec_id in specs:
            if spec_id in dereplicator_spec_to_compounds:
                clust_compounds |= dereplicator_spec_to_compounds[spec_id]
        if clust_compounds:
            cluster_to_compounds[cluster] = clust_compounds
            for spec_id in specs:
                if spec_id not in spec_putative_compounds:
                    spec_putative_compounds[spec_id] = clust_compounds
                else:
                    spec_putative_compounds[spec_id] |= clust_compounds

    for (clust1, clust2) in edges:
        if clust1 not in cluster_to_specs or clust2 not in cluster_to_specs:
            continue
        clust1_compounds = cluster_to_compounds[clust1] if clust1 in cluster_to_compounds else []
        clust2_compounds = cluster_to_compounds[clust2] if clust2 in cluster_to_compounds else []
        if clust2_compounds:
            for spec_id in cluster_to_specs[clust1]:
                if spec_id not in spec_putative_compounds:
                    spec_putative_compounds[spec_id] = clust2_compounds
                else:
                    spec_putative_compounds[spec_id] |= clust2_compounds
        if clust1_compounds:
            for spec_id in cluster_to_specs[clust2]:
                if spec_id not in spec_putative_compounds:
                    spec_putative_compounds[spec_id] = clust1_compounds
                else:
                    spec_putative_compounds[spec_id] |= clust1_compounds
    return spec_putative_compounds


def write_putative_compounds(cfg, spectra_fpaths, spec_putative_compounds):
    for idx, fpath in enumerate(spectra_fpaths):
        prepared_fpath = _get_fpath_by_idx(cfg, idx, splitext(fpath)[1])
        candidates_fpath = prepared_fpath + '.candidates'
        info('Writing candidate peptides for ' + fpath + ' to ' + candidates_fpath)
        with open(candidates_fpath, 'w') as f:
            spec_fname = basename(fpath)
            related_spec_ids = list(filter(lambda x: x[0] == spec_fname, spec_putative_compounds.keys()))
            for spec_id in related_spec_ids:
                scan = spec_id[1]
                unique_nlps = spec_putative_compounds[spec_id]
                num_nlps = len(unique_nlps)
                f.write(scan + ' ' + str(num_nlps))
                for nlp in unique_nlps:
                    f.write(' ' + nlp)
                f.write('\n')


def write_propagations(base_fname, PSMs, spec_to_cluster, cluster_to_spec, cluster_to_neighbours, edges,
                       spec_to_metadata, edge_to_metadata):
    # constants
    max_distance = 2
    #cluster_header_columns = '\t'.join(['ClusterId', 'ClusterSize', 'DeltaMZ', 'Cosine'])
    cluster_header_columns = '\t'.join(['ClusterId', 'ClusterSize'])
    indent = '\t'
    spectrum_header_columns = '\t'.join(['SpecFile', 'Scan', 'SpectrumMass', 'Retention'])

    def __print_level_header(stream, level):
        if level == 'root':  # special case:
            msg = 'Primary identification'
            additional_columns = PSM_ORF.header_line
            cur_indent = ''
        else:
            msg = 'Neighbours at distance %d edge(s)' % level
            additional_columns = spectrum_header_columns + '\n'
            cur_indent = indent
        stream.write(cur_indent + '### Level "{level}": {msg}\n'.format(**locals()))
        stream.write(cur_indent + '#' + cluster_header_columns + '\t' + additional_columns)

    def __print_cluster_info(stream, cluster_id, neighbour_cluster_id, verbose=True):
        cluster = cluster_to_spec[cluster_id]
        if neighbour_cluster_id is None:
            delta, cosine = 0.0, 1.0
        else:
            delta, cosine = edge_to_metadata[min(cluster_id, neighbour_cluster_id), max(cluster_id, neighbour_cluster_id)]
        for spec_id in cluster:
            pm, rt = spec_to_metadata[spec_id]
            stream.write(indent + '\t'.join(map(str, [cluster_id, len(cluster), #'%.2f' % delta, '%.2f' % cosine,
                                                      spec_id[0], spec_id[1], pm, rt])) + '\n')
            if not verbose:
                break

    def __get_psm_cluster(psm):
        spec_id = __get_spec_id(psm)
        cluster_id = spec_to_cluster[spec_id] if spec_id in spec_to_cluster else None
        return cluster_id

    def __print_psm_info(stream, psm, cluster_id):
        __print_level_header(stream, 'root')
        cluster = cluster_to_spec[cluster_id] if cluster_id in cluster_to_spec else []
        stream.write('\t'.join(map(str, [cluster_id, len(cluster), psm])) + '\n')

    def __breadth_first_search(stream, graph, root, verbose=True):
        if root is None:
            return []
        visited, queue = set(), collections.deque([(root, None, 0)])
        cur_distance = -1
        bfs_nodes = []
        while queue:
            vertex, prev_vertex, distance = queue.popleft()
            if distance > max_distance:
                break
            bfs_nodes.append(vertex)
            if distance > cur_distance:
                cur_distance = distance
                if verbose or cur_distance != 0:
                    __print_level_header(stream, cur_distance)
            if verbose or cur_distance != 0:
                __print_cluster_info(stream, vertex, prev_vertex, verbose=verbose)
            if vertex in graph:
                for neighbour in graph[vertex]:
                    if neighbour not in visited:
                        visited.add(neighbour)
                        queue.append((neighbour, vertex, distance + 1))
        return bfs_nodes

    def __draw_propagation(all_graphs, psm_nodes, root, title=''):
        def __contract_isotopes(graph, max_delta=2.5):
            mz_node_pairs = []
            for node in graph:
                mz, _ = spec_to_metadata[cluster_to_spec[node][0]]
                mz_node_pairs.append((mz, node))
            mz_node_pairs.sort()
            contracting_pair = mz_node_pairs[0]
            for cur_pair in mz_node_pairs[1:]:
                if cur_pair[0] - contracting_pair[0] <= max_delta:
                    if cur_pair[1] == root:  # special case -- don't remove root node
                        graph = nx.contracted_nodes(graph, cur_pair[1], contracting_pair[1])
                    else:
                        graph = nx.contracted_nodes(graph, contracting_pair[1], cur_pair[1])
                else:
                    contracting_pair = cur_pair
            return graph

        full_graph = nx.Graph()
        full_graph.add_edges_from(edges)
        psm_graph = full_graph.subgraph(psm_nodes)
        if len(psm_graph) > 1:
            psm_graph = __contract_isotopes(psm_graph)
        color_map = []
        labels = dict()
        for node in psm_graph:
            if node != root:
                color_map.append('green')
            else:
                color_map.append('blue')
            mz, rt = spec_to_metadata[cluster_to_spec[node][0]]
            labels[node] = 'mz: %.1f\nrt: %.1f' % (mz, rt)
        #graph_pos = nx.circular_layout(psm_graph)
        graph_pos = nx.spring_layout(psm_graph)
        nx.draw_networkx(psm_graph, graph_pos, node_color=color_map, with_labels=True,
                         node_size=1600, labels=labels, font_size=7)
        plt.title(title)
        plt.axis('off')
        all_graphs.savefig(plt.gcf())
        plt.close('all')

    def __make_plot_title(psm):
        spectrum_fname = basename(psm.spectrum_fpath)
        return 'Filename: {spectrum_fname}; Scan: {psm.scan_id} \n Compound: {psm.modified_seq}'.format(**locals())

    full_output_fpath = base_fname + '_detailed.txt'
    short_output_fpath = base_fname + '_summary.txt'
    graphs_output_fpath = base_fname + '.pdf'

    can_draw = True
    try:
        import networkx as nx
        import matplotlib
        matplotlib.use('Agg')  # non-GUI backend
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        all_graphs_pdf = PdfPages(graphs_output_fpath)
    except ImportError:
        can_draw = False
        warning('Cannot draw Spectral Network components: check that matplotlib and networkx Python libs are installed')

    num_skipped_psms = 0
    num_drawn_graphs = 0
    with open(full_output_fpath, 'w') as full_f:
        with open(short_output_fpath, 'w') as short_f:
            for psm in PSMs:
                cluster_id = __get_psm_cluster(psm)
                if cluster_id is None:
                    num_skipped_psms += 1
                    continue
                __print_psm_info(full_f, psm, cluster_id)
                __print_psm_info(short_f, psm, cluster_id)
                # TODO: think about running a single function with two streams (full and short)
                __breadth_first_search(full_f, cluster_to_neighbours, cluster_id)
                psm_nodes = __breadth_first_search(short_f, cluster_to_neighbours, cluster_id, verbose=False)
                if can_draw and len(psm_nodes) > 1:
                    num_drawn_graphs += 1
                    __draw_propagation(all_graphs_pdf, psm_nodes, cluster_id, __make_plot_title(psm))
    result_fpaths = [short_output_fpath, full_output_fpath]
    if can_draw:
        all_graphs_pdf.close()
        if num_drawn_graphs:
            result_fpaths.insert(0, graphs_output_fpath)
    if num_skipped_psms:
        warning('For %d out of %d significant PSMs, clusters in the Spectral Network were not found. '
                'Check that your Spectral Network is related to the same spectral dataset!' % (num_skipped_psms, len(PSMs)))
    return result_fpaths