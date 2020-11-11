#!/usr/bin/env python3
# Copyright (c) 2020 Oregon State University
# All Rights Reserved.
#
# Permission to use, copy, modify, and distribute this software and its
# documentation for educational, research and non-profit purposes, without
# fee, and without a written agreement is hereby granted, provided that
# the above copyright notice, this paragraph and the following three
# paragraphs appear in all copies.
#
# Permission to incorporate this software into commercial products may
# be obtained by contacting Oregon State University Office of Technology
# Transfer.
#
# This software program and documentation are copyrighted by Oregon State
# University. The software program and documentation are supplied "as is",
# without any accompanying services from Oregon State University. OSU does
# not warrant that the operation of the program will be uninterrupted or
# error-free. The end-user understands that the program was developed for
# research purposes and is advised not to rely exclusively on the program
# for any reason.
#
# IN NO EVENT SHALL OREGON STATE UNIVERSITY BE LIABLE TO ANY PARTY FOR
# DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
# INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS
# DOCUMENTATION, EVEN IF OREGON STATE UNIVERSITYHAS BEEN ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE. OREGON STATE UNIVERSITY SPECIFICALLY
# DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE AND
# ANY STATUTORY WARRANTY OF NON-INFRINGEMENT. THE SOFTWARE PROVIDED
# HEREUNDER IS ON AN "AS IS" BASIS, AND OREGON STATE UNIVERSITY HAS NO
# OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
# MODIFICATIONS.
from __future__ import print_function
import fileinput
import argparse
import os
import logging
import pprint as pp
import subprocess
import xml.etree.ElementTree as ElementTree
import re
import shutil
import gzip
import json
from math import ceil
from urllib.request import urlopen
from urllib.error import URLError
from rich.logging import RichHandler
from rich.console import Console
from rich.progress import (
    BarColumn,
    DownloadColumn,
    TextColumn,
    TransferSpeedColumn,
    TimeRemainingColumn,
    Progress,
    track,
)
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import as_completed
from functools import partial
from signal import signal, SIGPIPE, SIGINT, SIG_DFL
from .__version__ import __version__
signal(SIGPIPE, SIG_DFL)
signal(SIGINT, SIG_DFL)

NOT_WHITESPACE = re.compile(r'[^\s]')
OUTPUTS = ('fna', 'ffn', 'gff', 'gbk', 'faa')


class DownloadEdirect(argparse.Action):
    def __init__(self, nargs='?', **kw):
        super().__init__(nargs=nargs, **kw)

    def __call__(self, parser, namespace, values, option_string=None):
        init_logger(True)
        logger = logging.getLogger(__name__)
        if values:
            prefix = values
        else:
            prefix = '~/edirect'
        if '~' in prefix:
            prefix = os.path.expanduser(prefix)
        prefix = os.path.abspath(prefix)
        if os.path.basename(prefix) == 'edirect':
            dlpath = os.path.dirname(prefix)
        else:
            prefix = os.path.join(prefix, 'edirect')
            dlpath = prefix
        if not os.path.exists(dlpath):
            os.makedirs(dlpath)
        os.chdir(dlpath)
        out_file = os.path.join(dlpath, 'edirect.tar.gz')
        url = 'ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz'
        dl_file(url, out_file)
        logger.info(f'Installing edirect in {prefix}. Running setup.sh.')
        shutil.unpack_archive(out_file, dlpath, 'gztar')
        os.remove(out_file)
        subprocess.run([os.path.join(prefix, 'setup.sh')])
        logger.info('Add given directory {} to [green]$PATH[/green] to '
                    'continue.'.format(prefix), extra={'markup': True})
        logger.info('This script will find ~/edirect by default.')
        logger.info('Alternatively, set [green]$EDIRECT[/green] environment '
                    'variable.', extra={'markup': True})
        parser.exit()


def dl_file(url, out_file):
    logger = logging.getLogger(__name__)
    logger.debug('Downloading file {} to {}'.format(url, out_file))
    if not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
        with urlopen(url) as response, open(out_file, 'wb') as tarfh:
            shutil.copyfileobj(response, tarfh)


def dl_gzip(pbar, task_id, uri, filename, out, acc, swap):
    try:
        copy_url(pbar, task_id, uri, filename)
    except URLError:
        pbar.reset(task_id, start=False)
        uri = uri.replace('GCA_', 'GCF_')
        filename = filename.replace('GCA_', 'GCF_')
        copy_url(uri, task_id, uri, filename)
        # swap = True
    with gzip.open(filename, mode='rb') as gzfh,\
            open(out, 'wb') as outfh:
        shutil.copyfileobj(gzfh, outfh, 65536)
    # subprocess.run(['gunzip', dl_gz])
    # dl = dl_gz.replace('.gz', '')
    # os.rename(dl, out)
    if os.path.exists(out):
        # sys.stderr.write('\n')
        os.remove(filename)
    return out


def copy_url(pbar, task_id, uri, path) -> None:
    """Copy data from a url to a local file."""
    with urlopen(uri) as response:
        # This will break if the response doesn't contain content length
        pbar.update(task_id, total=int(response.info()["Content-length"]))
        with open(path, "wb") as dest_file:
            pbar.start_task(task_id)
            pbar.update(task_id, visible=True)
            for data in iter(partial(response.read, 32768), b""):
                dest_file.write(data)
                pbar.update(task_id, advance=len(data))
    pbar.update(task_id, visible=False)


def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x) and x != '-':
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


def edirect_dir(x):
    """
    'Type' for argparse - checks that edirect path is present.
    """
    init_logger(True)
    logger = logging.getLogger(__name__)
    if not os.path.exists(os.path.join(x, 'edirect.pl')):
        msg = ('Given path [bold magenta]{}[/bold magenta] does '
               '[bold red]NOT[/bold red] appear to be a valid edirect path.')
        logger.error(msg.format(x), extra={'markup': True})
        msg = ('Please specify a valid edirect path with '
               '[bold magenta]--edirect[/bold magenta] and try again.')
        logger.error(msg, extra={'markup': True})
        exit()
    return x


def init_logger(debug, logfile=''):
    logger = logging.getLogger()

    ch = RichHandler(rich_tracebacks=True)
    logger.setLevel(logging.DEBUG)
    if debug:
        ch.setLevel(logging.DEBUG)
    else:
        ch.setLevel(logging.INFO)
    logger.addHandler(ch)

    if logfile:
        fh = RichHandler(console=Console(file=open(logfile, 'a')))
        logger.addHandler(fh)


def run_argparse():
    parser = argparse.ArgumentParser(
        description='Downloads assemblies & annotations from NCBI.')
    subparsers = parser.add_subparsers(
        dest='intype', help='Choose from this list of input types.'
    )
    parser.add_argument(
        '--debug', help='Turn on debugging messages.', action='store_true'
    )
    parser.add_argument(
        '--version', action='version', version=__version__
    )
    parser.add_argument(
        '--dledirect', help='Download edirect to given location. [~/edirect]',
        action=DownloadEdirect, type=str, nargs='?'
    )

    # Add input types
    organism = subparsers.add_parser('organism',
                                     help='Valid NCBI organism or taxids.')
    assembly_ids = subparsers.add_parser('assembly_ids',
                                         help='Valid NCBI assembly IDs.')
    nuccore_ids = subparsers.add_parser('nuccore_ids',
                                        help='Valid NCBI nucleotide '
                                             'accessions.')
    json_input = subparsers.add_parser('json_input',
                                       help='Valid NCBI JSON docsums.')

    p_assemlinks = [organism, assembly_ids, nuccore_ids]
    all_p = [organism, assembly_ids, nuccore_ids, json_input]

    organism.add_argument(
        'query', help='Valid NCBI organism text term or ID', type=str
    )
    organism.add_argument(
        '--type', help='Input is text term (default) or ID', type=str,
        choices=['text', 'ID'], default='text'
    )

    assembly_ids.add_argument(
        'infile', help='Input file with NCBI assembly IDs; "-" for stdin',
        type=extant_file
    )
    assembly_ids.add_argument(
        '--type', help='Input is Accession (default) or ID', type=str,
        choices=['acc', 'uid'], default='acc'
    )

    nuccore_ids.add_argument(
        'infile', help='Input file with NCBI nuccore IDs; "-" for stdin',
        type=extant_file
    )

    json_input.add_argument(
        'jsonfile', help='Input JSON file with docsums; "-" for stdin',
        type=extant_file, nargs='+'
    )
    json_input.add_argument(
        '--function', help='Download metadata and/or genomes. [metadata]',
        type=str, default=['metadata'], choices=['metadata', 'genomes'],
        nargs='+'
    )

    for p in p_assemlinks:
        p.add_argument(
            '--function', help='check counts, download metadata,'
                               ' or genomes. [check]', type=str,
            default=['check'], choices=['check', 'metadata', 'genomes'],
            nargs='+'
        )

    for p in all_p:
        p.add_argument(
            '--annotation', help='Require annotation? False by default,'
                                 ' True if gbk/faa/ffn requested',
            action='store_true', default=False
        )
        p.add_argument(
            '--metadata_append', help='Append to metadata, not overwrite.',
            action='store_true', default=False
        )
        p.add_argument(
            '--typestrain', help='Only download type strains.',
            action='store_true', default=False
        )
        p.add_argument(
            '--keepmulti', help='By default, genomes from large multi-isolate'
            'studies are removed.', action='store_true', default=False
        )
        p.add_argument(
            '--force', help='Force download attempt of low-quality genomes.',
            action='store_true', default=False
        )
        p.add_argument(
            '-f', '--outformat', help='Output file prefix. [full]', type=str,
            default='full', choices=['abbr', 'full', 'strain']
        )
        p.add_argument(
            '-o', help='Output file types.', type=str,
            choices=['fna', 'ffn', 'gff', 'gbk', 'faa', 'all'], nargs='+'
        )
        p.add_argument(
            '--edirect', help='Path to edirect directory.', type=edirect_dir
        )
        p.add_argument(
            '--debug', help='Turn on debugging messages.', action='store_true'
        )

    args = parser.parse_args()
    init_logger(args.debug, logfile='get_assemblies.log')
    if not args.intype:
        parser.print_help()
        exit()
    # if args.debug:
    #     pp.pprint(args)
    return args


def decode_stacked_json(document, pos=0, decoder=json.JSONDecoder()):
    while True:
        match = NOT_WHITESPACE.search(document, pos)
        if not match:
            return
        pos = match.start()

        try:
            obj, pos = decoder.raw_decode(document, pos)
        except json.JSONDecodeError:
            # do something sensible if there's some error
            raise
        yield obj


def eutil_critical(name):
    logger = logging.getLogger(__name__)
    logger.critical(
        f'Unable to find eutil: {name}'
    )
    logger.critical(
        'Use --dledirect to download the edirect files and try again.'
    )
    exit(1)


def validate_inputs(args):
    logger = logging.getLogger(__name__)

    if 'genomes' in args.function and 'metadata' not in args.function:
        if args.intype != 'json_file':
            args.function.append('metadata')

    if 'genomes' in args.function and not args.o:
        logger.critical(
            '--function genomes set, but no outputs chosen using --outputs'
        )
        logger.critical(
            'Check settings (add at least one output) and try again.'
        )
        exit(1)
    else:
        if args.o:
            if 'all' in args.o:
                args.o = OUTPUTS
            else:
                args.o = set(args.o)

    if args.o:
        for opt in args.o:
            if opt in ('ffn', 'gbk', 'faa'):
                args.annotation = True
                break

    user_edir = os.path.expanduser(os.path.join('~', 'edirect'))
    if args.edirect:
        edirect_dir = args.edirect
    elif 'EDIRECT' in os.environ and os.path.exists(os.environ['EDIRECT']):
        edirect_dir = os.environ['EDIRECT']
    elif os.path.exists(user_edir):
        edirect_dir = user_edir
    else:
        logger.debug('No edirect dir found.')
        edirect_dir = ''

    utils = ['esearch', 'efetch', 'xtract', 'efilter', 'elink', 'epost']
    exes = {}
    for util in utils:
        if shutil.which(util):
            exes[util] = util
        elif os.path.exists(edirect_dir):
            check = os.path.join(edirect_dir, util)
            if os.path.isfile(check):
                exes[util] = check
        else:
            eutil_critical(util)

    # logger.debug(pp.pformat(args))
    logger.debug(pp.pformat(args))
    logger.debug(pp.pformat(exes))
    return(args, exes)


def search_query(esearch, efilter, qtype, query):
    logger = logging.getLogger(__name__)
    # returns list of assembly IDs
    # input is organism search or taxid
    command = []
    if qtype == 'text':
        query = f'{query}[ORGN]'
    elif qtype == 'ID':
        query = f'txid{query}[ORGN]'

    query = query + ' AND latest [PROP] NOT '\
        'partial-genome-representation [PROP]'

    command = [esearch, '-db', 'assembly', '-query', query]

    # logger.debug('Running command:\n' + ' '.join(command))
    logger.debug('Running command:\n' + ' '.join(command))
    try:
        search_res = subprocess.run(command, stdout=subprocess.PIPE, text=True,
                                    check=True)
    except subprocess.CalledProcessError as e:
        logger.critical('There was an error running esearch.')
        logger.critical(f'Check the error: "{e}" and try again.')
        exit(1)

    return(search_res.stdout)


def get_assem_links(epost, infile, qtype):
    # returns list of assembly IDs
    # input is list of assem ids
    logger = logging.getLogger(__name__)
    assem_ids = [x.strip() for x in fileinput.input(infile) if x]
    logger.debug(pp.pformat(assem_ids))

    command = [epost,
               '-db', 'assembly',
               '-format', f'{qtype}']

    # command = ' '.join(command)
    try:
        output = subprocess.run(command, stdout=subprocess.PIPE, shell=False,
                                text=True, input='\n'.join(assem_ids),
                                check=True)
    except subprocess.CalledProcessError as e:
        logger.critical('There was an error running epost.')
        logger.critical(f'Check the error: "{e}" and try again.')
        exit(1)
    return(output.stdout)


def convert_nuc_to_assem(elink, infile):
    # returns list of assembly IDs
    # input is a list of nucleotide IDs
    logger = logging.getLogger(__name__)
    nuccore_ids = [x.strip() for x in fileinput.input(infile) if x]
    logger.debug(pp.pformat(nuccore_ids))

    command = [elink,
               '-target', 'assembly',
               '-format', 'acc',
               '-db', 'nuccore',
               '-id',
               ','.join(nuccore_ids)]

    # command = ' '.join(command)
    try:
        output = subprocess.run(command, stdout=subprocess.PIPE, shell=False,
                                text=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.critical('There was an error running elink.')
        logger.critical(f'Check the error: "{e}" and try again.')
        exit(1)

    return(output.stdout)


def check_count(assem_links, query='stdin'):
    logger = logging.getLogger(__name__)
    try:
        entrez_direct = ElementTree.fromstring(assem_links)
    except ElementTree.ParseError:
        logger.critical(
            'No output from edirect search.'
        )
        logger.critical(
            'Check input(s) and try again.'
        )
    try:
        count = entrez_direct.find('Count').text
    except AttributeError:
        with open('no_hits.txt', 'a') as nohitfh:
            nohitfh.write(f'{query}\n')
        logger.critical(
            f'No genomes found to download with query {query}.'
        )
        logger.critical(
            'Check your query/taxid and try again.'
        )
        logger.critical(
            'Examine no_hits.txt for a record of genomes with no hits.'
        )
        exit(1)
    else:
        logger.info(
            f'Found {count} genome(s) to download.'
        )
        logger.info(
            f'Expect {int(count)*5} MB to {int(count)*7} MB of data pending '
            'the chosen file types for download.'
        )
    # return(count)


def chunks(uidl, n):
    """Yield successive n-sized chunks from list."""
    for i in range(0, len(uidl), n):
        yield uidl[i:i + n]


def fetch_docsums(efetch, assem_links):
    logger = logging.getLogger(__name__)
    command = [f'{efetch}', '-format', 'uid']
    try:
        uids = subprocess.run(' '.join(command), stdout=subprocess.PIPE,
                              shell=True, text=True, input=assem_links,
                              check=True)
    except subprocess.CalledProcessError as e:
        logger.critical('There was an error running efetch.')
        logger.critical(f'Check the error: "{e}" and try again.')
        exit(1)
    # logger.debug('Found these uids:\n' +
    #                   uids.stdout)
    uid_list = sorted(uids.stdout.splitlines())

    outputs = []
    i = 0
    nchunk = ceil(len(uid_list)/500)
    if nchunk > 0:
        units = 'secs'
        length = nchunk*10
        if length > 60:
            units = 'mins'
            length = length/60
            if length > 60:
                units = 'hours'
                length = length/60
        logger.info('With {} chunks (500 ids per), this will take around '
                    '{} {}.'.format(nchunk, int(length), units))
    for chunk in track(chunks(uid_list, 500), 'chunk', nchunk):
        command = [f'{efetch}',
                   '-format', 'docsum',
                   '-mode', 'json',
                   '-db', 'assembly',
                   '-id', ','.join(chunk)]
        # command += ','.join(chunk)
        # logger.debug('Running this command:\n' + ' '.join(command))
        logger.debug('Running this command: ' + ' '.join(command))
        jsondata = []

        for j in range(5):
            try:
                output = subprocess.run(command, stdout=subprocess.PIPE,
                                        shell=False, text=True, check=True)
            except subprocess.CalledProcessError as e:
                logger.critical('There was an error running efetch.')
                logger.critical(f'Check the error: "{e}" and try again.')
                exit(1)

            try:
                jsondata = json.loads(output.stdout)
            except json.JSONDecodeError:
                msg = 'No JSON output detected. Retrying. ({})'
                logger.warning(msg.format(j+1))
            else:
                break
        if not jsondata:
            msg = 'Unable to download chunk from NCBI.'
            logger.critical(msg)
            msg = 'Check settings and try again later.'
            logger.critical(msg)
            exit(-1)
        with open(f'docsums{i}.json', 'w') as docsumfh:
            json.dump(jsondata, docsumfh, indent=4)
            docsumfh.write('\n')
        outputs.append(jsondata)
        i += 1

    # with open('docsums.json', 'w') as fh:
    #     fh.write(outputs)
    #     fh.write('\n')

    return(outputs)


def get_json_from_file(jsonfiles):
    # returns list of assembly IDs
    # input is list of assem ids
    docsums = []
    for jsonfile in jsonfiles:
        with open(jsonfile, 'r') as fh:
            docsum = json.load(fh)
        docsums.append(docsum)
    return(docsums)


def json_error(datatype, infotype, argument):
    # Print can't find datatype in argument
    logger = logging.getLogger(__name__)
    logger.critical(f'I cannot find the {datatype} from this {infotype}.')
    logger.critical(f'{infotype} -- {argument}')
    logger.critical(
        f'Check that there are {datatype}s returned and try again')
    exit(-1)


def debugging(datatype, infotype, argument):
    # Print can't find datatype in argument
    logger = logging.getLogger(__name__)
    logger.debug(f'I cannot find the {datatype} from this {infotype}.')
    logger.debug(f'{infotype} -- {argument}')


def remove_invalid_characters(string):
    string = string.replace('[', '')
    string = string.replace(']', '')
    string = string.replace('/', '_')
    string = string.replace('(', '')
    string = string.replace(')', '')
    string = string.replace(':', '_')
    string = string.replace(';', '')
    string = string.replace('%', '')
    return string


def replace_spaces(string):
    string = string.replace(' ', '_')
    return string


def get_prefix(outformat, name, strain, assem_name):
    # Get valid genus and species
    # Omit Candidatus, if present
    logger = logging.getLogger(__name__)
    name = remove_invalid_characters(name)
    strain = remove_invalid_characters(strain)

    name_list = name.split(' ')
    if name_list[0] == 'Candidatus' or name_list[0] == 'uncultured':
        if len(name_list) > 2:
            genus = name_list[1]
            species = name_list[2]
            strain_idx = 3
        else:
            genus = name_list[1]
            species = 'sp'
            strain_idx = 2
    else:
        try:
            genus = name_list[0]
        except IndexError:
            genus = ''
            logger.error(f'Unable to find genus for {name}')
        try:
            species = name_list[1]
        except IndexError:
            species = ''
            logger.error(f'Unable to find species for {name} {strain}')
        strain_idx = 2

    if genus in strain and 'sp.' not in name:
        warn = True
        logger.debug(f'Seemingly malformed organism name {name} {strain}')
        logger.debug(f'Strain before {strain}')
        try:
            strain = strain.split(' ',  2)[2]
        except IndexError:
            strain = '-'
        logger.debug(f'Strain after {strain}')
    else:
        warn = False

    if strain == '-':
        if len(name_list) == strain_idx + 1:
            strain = name_list[strain_idx]
        else:
            strain = ' '.join(name_list[strain_idx:])

    # Remove spaces
    genus = replace_spaces(genus)
    species = replace_spaces(species)
    strain = replace_spaces(strain)

    if not strain:
        strain = assem_name

    if 'substr.' in strain:
        strain = strain.replace('substr.', 'substr')

    if outformat == 'strain' or not species:
        prefix = strain
    elif outformat == 'full' or 'sp.' in species:
        species = species.replace('sp.', 'sp')
        prefix = '_'.join([genus, species, strain])
    elif outformat == 'abbr':
        prefix = '_'.join([genus[0], species, strain])

    if warn:
        logger.warning(f'Check out this prefix for accuracy: {prefix}')

    return prefix


def extract_metadata(force, metadata_append, outformat, typestrain, annotation,
                     docsums, keepmulti):
    # Return dl_mapping dict at the end
    logger = logging.getLogger(__name__)
    dl_mapping = {}
    # d = json.loads(docsums)

    json_keys = ['uid', 'assemblyname', 'speciesname', 'wgs',
                 'submitterorganization', 'assemblystatus',
                 'fromtype', 'lastupdatedate', 'coverage', 'contign50']
    special_keys = ['isolate/strain', 'sequence_type',
                    'accession_type', 'accession']
    prefixes = {}
    metadata = []
    metadata_file = 'metadata.tab'
    if metadata_append:
        open_mode = 'a'
        if not os.path.exists(metadata_file):
            metadata.append(json_keys + special_keys + ['prefix'])
    else:
        open_mode = 'w'
        metadata.append(json_keys + special_keys + ['prefix'])

    # for d in decode_stacked_json(docsums):
    for d in track(docsums, 'docsums', len(docsums)):
        try:
            json_data = d['result']
        except KeyError:
            json_error('result', 'json_docsums', docsums)

        try:
            uids = json_data['uids']
        except KeyError:
            json_error('uid', 'json_docsums', docsums)

        # Keep track of strain names so we don't have dups

        for uid in uids:
            # annotation_type = 'genbank'
            annotation_type = 'refseq'
            line = []
            skip = False
            # Standard keys; only str/int return
            for json_key in json_keys:
                try:
                    item = json_data[uid][json_key]
                except KeyError:
                    debugging(json_key, 'uid', uid)
                    line.append('-')
                else:
                    if item:
                        line.append(str(item))
                    else:
                        line.append('-')
                    if json_key == 'fromtype':
                        if not item and typestrain:
                            skip = True
                    if json_key == 'speciesname':
                        name = item
                    elif json_key == 'assemblyname':
                        assem_name = item.replace(' ', '_')
                    elif json_key == 'assemblystatus':
                        if item == 'Complete Genome':
                            annotation_type = 'genbank'
            if skip:
                logger.debug(
                    f'{uid} - {name} {assem_name} not from type and '
                    '--typestrain set. Skipping...'
                )
                continue

            # Special values; need more digging through list/dict
            # sequence_type
            sequence_type = ''
            try:
                exclfromrefseq = json_data[uid]['exclfromrefseq']
            except KeyError:
                debugging('exclfromrefseq', 'uid', uid)
                exclfromrefseq = []
            skip = False
            reason = ''
            # https://www.ncbi.nlm.nih.gov/assembly/help/anomnotrefseq/
            for item in exclfromrefseq:
                annotation_type = 'genbank'
                if item == 'derived from single cell':
                    sequence_type = 'SAG'
                if item == 'derived from metagenome':
                    sequence_type = 'metagenome'
                if item == 'derived from environmental source':
                    sequence_type = 'environmental'
                if item == 'from large multi-isolate project':
                    sequence_type = 'multi-isolate'
                    if not keepmulti:
                        skip = True
                        reason = item
                if item in ('low contig N50', 'many frameshifted proteins',
                            'low quality sequence', 'genome length too large',
                            'contaminated', 'abnormal gene to sequence ratio',
                            'fragmented assembly', 'genome length too small',
                            'low gene count', 'missing ribosomal protein '
                            'genes', 'missing rRNA genes', 'missing tRNA '
                            'genes', 'partial', 'untrustworthy as type',
                            'chimeric', 'hybrid', 'misassembled', 'sequence '
                            'duplications', 'unverified source organism'):
                    skip = True
                    reason = item

            # Strain
            strain = ''
            isolate = ''
            sub_value = ''
            try:
                biosource = json_data[uid]['biosource']
            except KeyError:
                strain = '-'
                debugging('biosource', 'uid', uid)
            else:
                try:
                    isolate = biosource['isolate']
                except KeyError:
                    isolate = ''
                    debugging('isolate', 'biosource', uid)
                for item in biosource['infraspecieslist']:
                    try:
                        sub_type = item['sub_type']
                    except KeyError:
                        continue
                    else:
                        if sub_type == 'strain':
                            try:
                                sub_value = item['sub_value']
                            except KeyError:
                                debugging('strain', 'biosource', uid)
                                sub_value = ''
                if sub_value:
                    strain = sub_value
                elif isolate:
                    strain = isolate
                else:
                    debugging('strain_name', 'biosource', uid)
                    strain = '-'
            # try:
            #     seen[strain] += 1
            # except KeyError:
            #     seen[strain] = 0
            # else:
            #     strain = f'{strain}_{seen[strain]}'
            #     logger.debug(f'Adding number to strain {strain}.')
            line.append(strain)

            # sequence_type continued
            if skip:
                if force:
                    logger.warning(
                        f'Assembly {uid} - {strain} is being skipped due '
                        f'to "{reason}". --force enabled; attempting download.'
                    )
                else:
                    logger.warning(
                        f'Assembly {uid} - {strain} is being skipped due '
                        f'to "{reason}".'
                    )
                    continue

            try:
                propertylist = json_data[uid]['propertylist']
            except KeyError:
                annotation_type = False
                debugging('propertylist', 'uid', uid)
            else:
                if 'latest' not in propertylist:
                    if not force:
                        logger.debug(
                            f'{strain} -- {uid} not latest. Skipping.'
                        )
                        continue
                if 'partial-genome-representation' in propertylist:
                    if not force:
                        logger.debug(
                            f'{strain} -- {uid} is partial genome. Skipping.'
                        )
                        continue

            if annotation:
                logger.debug(
                    f'{strain} -- annotation type: {annotation_type}'
                )
                if annotation_type == 'genbank':
                    if 'genbank_has_annotation' in propertylist:
                        pass
                    elif 'refseq_has_annotation' in propertylist:
                        logger.debug(
                            f'Switching to refseq annotation for {strain}'
                        )
                        annotation_type = 'refseq'
                    else:
                        logger.critical(
                            f'Unable to find annotation for {strain}.'
                            ' Skipping...'
                        )
                        continue
                elif annotation_type == 'refseq':
                    if 'refseq_has_annotation' in propertylist:
                        pass
                    elif 'genbank_has_annotation' in propertylist:
                        logger.debug(
                            f'Switching to genbank annotation for {strain}'
                        )
                        annotation_type = 'genbank'
                    else:
                        logger.critical(
                            f'Unable to find annotation for {strain}.'
                            ' Skipping...'
                        )
                        continue

            if not sequence_type:
                sequence_type = 'cell culture'
            line.append(sequence_type)

            # Accession
            accession = ''
            acc_type = 'ftppath_'
            try:
                synonym = json_data[uid]['synonym']
            except KeyError:
                logger.warning(f'Unable to find accession for {uid}.')
                logger.warning(f'{uid} will be skipped.')
                continue
            try:
                accession = synonym[annotation_type]
            except KeyError:
                accession = ''
            else:
                acc_type += annotation_type
            if not accession:
                logger.warning(f'Unable to find accession for {uid}. '
                               'Skipping...')
                continue
            else:
                line.append(annotation_type)
                line.append(accession)

            # Get file prefix
            prefix = get_prefix(outformat, name, strain, assem_name)
            if prefix in prefixes:
                logger.warning(f'Prefix "{prefix}" has already been seen, '
                               f'adding assembly name "{assem_name}" to end.')
                logger.debug('Already seen assembly ID: "{}", new ID "{}"'
                             .format(prefixes[prefix], accession))
                prefix = '_'.join([prefix, assem_name])
            else:
                prefixes[prefix] = accession
            line.append(prefix)

            # Get FTP path
            ftp_path = ''
            try:
                ftp_path = json_data[uid][acc_type]
            except KeyError:
                logger.warning(f'Unable to find ftppath for {accession}.')
                logger.warning(f'{accession} will be skipped.')
                continue

            dl_mapping[accession] = {'ftp_path': ftp_path,
                                     'prefix': prefix,
                                     'assem_name': assem_name,
                                     'swap': False}
            metadata.append(line)
    with open(metadata_file, open_mode) as metadatafh:
        for line in metadata:
            metadatafh.write('\t'.join(line) + '\n')

    return dl_mapping


def download_genomes(o, dl_mapping):
    logger = logging.getLogger(__name__)
    filemap = {'faa': 'protein.faa.gz',
               'fna': 'genomic.fna.gz',
               'gbk': 'genomic.gbff.gz',
               'gff': 'genomic.gff.gz',
               'ffn': 'cds_from_genomic.fna.gz'}

    nacc = len(dl_mapping)
    logger.info('Downloading {} file(s).'.format(nacc * len(o)))

    # Generate urls to download
    urls = []
    for acc in dl_mapping:
        dl_base = '_'.join([acc, dl_mapping[acc]['assem_name']])
        dl_base = dl_base.replace(',', '')
        for ft in o:
            # Get output filename
            out = dl_mapping[acc]['prefix'] + '.' + ft
            if not os.path.exists(out):
                uri = dl_mapping[acc]['ftp_path'] + \
                    '/' + dl_base + '_' + filemap[ft]
                if dl_mapping[acc]['swap']:
                    uri = uri.replace('GCA_', 'GCF_')
                logger.debug(f'URL {uri} for {acc} generated.')
                urls.append((uri, out, acc, dl_mapping[acc]['swap']))
            else:
                logger.info(
                    f'{out} already found. Skipping.'
                )

    progress = Progress(
        TextColumn("[bold blue]{task.fields[filename]}", justify="right"),
        BarColumn(bar_width=None),
        "[progress.percentage]{task.percentage:>3.1f}%",
        "•",
        DownloadColumn(),
        "•",
        TransferSpeedColumn(),
        "•",
        TimeRemainingColumn(),
    )

    results = []
    with progress as pbar:
        with ThreadPoolExecutor(max_workers=4) as pool:
            for uri, out, acc, swap in urls:
                logger.debug(
                    f'Downloading {acc} to {out}.'
                )
                filename = uri.split("/")[-1]
                task_id = pbar.add_task('download', filename=out, start=False,
                                        visible=False)
                future = pool.submit(
                    dl_gzip, progress, task_id, uri, filename, out, acc, swap)
                results.append(future)
                # dl_mapping[acc]['swap'] = dl_gzip(
                #     pbar, task_id, uri, filename, out, acc, swap)

    for future in as_completed(results):
        output = future.result()
        if os.path.exists(output):
            logger.info(f'{output} successfully downloaded.')
        else:
            logger.error(f'{output} was unable to be downloaded.')


def main():
    args = run_argparse()
    args, exes = validate_inputs(args)

    # assem_ids is a list of NCBI assembly IDs
    if args.intype == 'organism':
        assem_links = search_query(exes['esearch'], exes['efilter'],
                                   args.type, args.query)
    elif args.intype == 'assembly_ids':
        assem_links = get_assem_links(exes['epost'], args.infile, args.type)
    elif args.intype == 'nuccore_ids':
        assem_links = convert_nuc_to_assem(exes['elink'], args.infile)
    try:
        with open('.assemlinks', 'wt') as afh:
            afh.write(assem_links)
    except UnboundLocalError:
        pass

    if args.intype in ['organism',
                       'assembly_ids',
                       'nuccore_ids']:
        try:
            check_count(assem_links, args.query)
        except AttributeError:
            check_count(assem_links)

    if 'metadata' in args.function:
        if args.intype == 'json_input':
            docsums = get_json_from_file(args.jsonfile)
        else:
            docsums = fetch_docsums(exes['efetch'], assem_links)
        dl_mapping = extract_metadata(args.force, args.metadata_append,
                                      args.outformat, args.typestrain,
                                      args.annotation, docsums, args.keepmulti)

    if 'genomes' in args.function:
        download_genomes(args.o, dl_mapping)


if __name__ == '__main__':
    main()
