VALUE_MAX_SCAN = 2000000
VALUE_APP_SCAN = 100000
VALUE_ILLEGAL = -716
VALUE_MAX_VALUE = 2000000

MS2_OUTPUT_TYPE = {
    'RAW': 0,
    'MS2_PROCESSED': 1
}

PSM_TYPE = {
    'RAW_PSM': 'rpsm',
    'FINE_SCORING_PSM':'fpsm',
    'SVM_PSM':'cxsvmpsm',
    'SITE_SCORING_PSM':'spsm',
    'SITE_FDR_PSM':'sfdrpsm',
}

PSM_OUTPUT_TYPE = {
    'ALL': 0,
    'TOP1_UNFILTERED': 1,
    'TOP1_FILTERED': 2,
}

GLOBAL_DIGEST_NAME = {
    'SPECIFIC': 0,
    'SEMI_SPECIFIC': 1,
    'UNSPECIFIC': 2
}

TOLERANCE_TYPE = {
    'ABSOLUTE_TOLERANCE': 0,
    'RELATIVE_TOLERANCE': 1
}

EXPIRATION_TIME = {'Year': 2419, 'Month': 12, 'Day': 31}

INFO_TO_USER_Staff = (
    '\n[ComMap] Copyright \u00A9 2021. All rights reserved.',
    '\n[ComMap] ComMap is expired! Please send e-mail for the new version.',
    '\n[ComMap] Warning! The current license will expired in 7 days. Please send e-mail for the new version.',
    '\n[ComMap] Writing config file in the folder...',
    '\n[ComMap] Finished!',
    '\n[ComMap] Writing config file...',
    '\n[ComMap] Run ComMap with this command: python ComMap [parameter file]',
    '\n[ComMap] ComMap whole workflow.',)

CFLOW6_INFORMATION = (
    #0
    '[ComMap] Start to read fasta file and cross-link result file.',
    #1
    '[ComMap] Start to download Uniprot information.',
    #2
    '[ComMap] End downloading Uniprot information.',
    #3
    '[ComMap] Start to download Protein Data Bank information.',
    #4
    '[ComMap] End downloading Protein Data Bank information.',
    #5
    '[ComMap] Start sequence alignment.',
    #6
    '[ComMap] End sequence alignment.',
    #7
    '[ComMap] Start structure distance calculation.',
    #8
    '[ComMap] End structure distance calculation.',
    #9
    '[ComMap] Organize and output.',
    #10
    '[ComMap] Preparing local PDB structures.',
    #11
    '[ComMap] Finish preparing local PDB structures.',
    #12
    '[ComMap] ERROR NO PDB STRUCTURE FOUNDED.',
    #13
    '[ComMap] Start quality control.',
    #14
    '[ComMap] ERROR READ MAP FILE FAILED.',
    #15
    '[ComMap] .',
)

# HTTP & FTP response result,Don't change it!
RESPONSE_RESULT = {
    'Fail': 0,
    'Success': 1
}

# The base URL of Uniprot database, Don't change it!
UNIPROT_BASE_URL = 'http://www.uniprot.org/uniprot/'

# The return file type from Uniprot database, Don't change it!
UNIPROT_FILE_TYPE = '.xml.gz'

# The base URL of Protein Data Bank database, Don't change it!
PDBJ_PDB_FORMAT_BASE_URL = 'ftp://ftp.pdbj.org/pub/pdb/data/structures/divided/pdb/'
PDBJ_XML_FORMAT_BASE_URL = 'ftp://ftp.pdbj.org/pub/pdb/data/structures/all/XML/'
PDBJ_CIF_FORMAT_BASE_URL = 'ftp://ftp.pdbj.org/pub/pdb/data/structures/divided/mmCIF/'

# The return file type from Protein Data Bank database, Don't change it!
PDBJ_PDB_FILE_TYPE = '.ent.gz'
PDBJ_XML_FILE_TYPE='.xml.gz'
PDBJ_CIF_FILE_TYPE='.cif.gz'

MAX_CONCURRENT_NUMBER=20

AMINO_ACIDS_DIC = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V',
    'UNK': 'X',
    'SEC': 'U',
}