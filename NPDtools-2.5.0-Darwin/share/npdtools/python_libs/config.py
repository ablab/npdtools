import platform

# set 'False' for Release which will hide very advanced options and defaults
# Note: GNPS-release requires DEVELOPER_VERSION == True
DEVELOPER_VERSION = False

MASS_PROTON = 1.00728

CRITICAL_SPECTRA_SCAN_LEN = 5000
DEFAULT_THREADS = 4  # Maximum default value, it may be decreased in some cases inside Dereplicator/RiPPquest wrappers
DEFAULT_P_VALUE_THRESHOLD = 1e-10
ACCURATE_P_VALUE_THRESHOLD = 1e-10

LOGGER_METAMINER_NAME = 'metaminer'
LOGGER_DEREPLICATOR_NAME = 'dereplicator'

DEFAULT_DATABASE = 'combined'

if platform.system() == 'Darwin':
    platform_name = 'osx'
else:
    platform_name = 'linux'
    # if platform.architecture()[0] == '64bit':
    #     platform_name = 'linux_64'
    # else:
    #     platform_name = 'linux_32'

dereplicator_p_value_threshold = DEFAULT_P_VALUE_THRESHOLD
rippquest_p_value_threshold = DEFAULT_P_VALUE_THRESHOLD

global_config = None

sequence_extensions_dict = {'fasta': ['.fna', '.fasta', '.fa'], 'antismash': ['.gbk'], 'boa': ['.txt']}
allowed_sequence_extensions = [item for sublist in sequence_extensions_dict.values() for item in sublist]
allowed_spectrum_extensions = ['.mzML', '.mzXML', '.mzdata', '.mz5', '.mgf', '.ms1', '.ms2', '.cms1', '.cms2', '.raw']

# RiPPquest Product Classes
prod_classes = ['formylated', 'glycocin', 'lantibiotic', 'lap', 'lassopeptide', 'linaridin', 'proteusin',
                'cyanobactin', 'methanobactin']
excluded_prod_class = ['formylated', 'lantibiotic_heavy']  # classes which are excluded from 'all product classes' runs

# Dereplicator+ constants
DEFAULT_DER_PLUS_SIG_MIN_SCORE = 12
DER_PLUS_ALL_MIN_SCORE_RATE = 3. / 4
der_plus_score_threshold = DEFAULT_DER_PLUS_SIG_MIN_SCORE
USE_SCORE_ONLY = False

# best spectra params for --extract-best-spectra (p-value < BS_MAX_P_VALUE, fdr < BS_MAX_FDR)
# Note: BS_MAX_P_VALUE should be <= P_VALUE_THRESHOLD
BS_MAX_P_VALUE = 1e-10
BS_MAX_FDR = 0.01

# RiPPquest antiSMASH processing
neighbours_distance = 10000  # for looking around secondary metabolite and determining its Production Class

# for copying in work directory
config_dirs = ['Fragmentation_rule', 'configs']
