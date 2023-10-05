#!/usr/bin/env python2.7

import argparse
import os
import tarfile
import shlex
import subprocess
import logging



class ExtLog:

    def __init__(self, name, path, level):
        self.__name__ = name
        self.__path__ = path
        if not os.path.exists(self.__path__):
            os.makedirs(self.__path__)
        self.log = logging.getLogger(name)
        self.log_file = os.path.join(self.__path__, '%s.log' % self.__name__)
        formatter = logging.Formatter(
            '%(levelname)s : %(asctime)s : %(message)s')
        fileHandler = logging.FileHandler(self.log_file, mode='w')
        fileHandler.setFormatter(formatter)
        streamHandler = logging.StreamHandler()
        streamHandler.setFormatter(formatter)
        self.log.setLevel(level)
        self.log.addHandler(fileHandler)
        self.log.addHandler(streamHandler)


def mkdir_sv_data(tool, path, log):
    log.log.info('Create %s for tool %s' % (path, tool))
    if not os.path.exists(path):
        os.makedirs(path)
    else:
        log.log.warning('Folder %s (for %s) already exists!' % (path, tool))


def symlink_sv_data(tool, sv_file, path, log):
    sv_path, sv_file_name = os.path.split(sv_file)
    sv_link = os.path.join(path, sv_file_name)
    log.log.info('Tool %s, symlink %s to %s' % (tool, sv_file, sv_link))
    os.symlink(sv_file, sv_link)
    return sv_link


def main():
    parser = argparse.ArgumentParser(description='Merge various sv caller')
    parser.add_argument('--run-id',  dest='run_id',
                        help='Sample id, identifier of the run', required=True)
    parser.add_argument('--delly',  dest='delly',
                        help='DELLY output VCF',  required=True)
    parser.add_argument('--snowman',  dest='snowman',
                        help='SNOWMAN output VCF',  required=True)
    parser.add_argument('--dranger',  dest='dranger',
                        help='DRANGER output VCF',  required=True)
    parser.add_argument('--brass',  dest='brass',
                        help='BRASS output VCF',
                        required=True)

    args = parser.parse_args()
    current_dir = os.getcwd()
    output_dir = os.path.expanduser('~')
    #log_dir = os.path.join(current_dir, 'logs')
    log_dir = os.path.join(output_dir, 'logs')
    log = ExtLog('run_dockstore', log_dir, level=logging.INFO)
    log.log.info('Output results in folder %s' % output_dir)
    input_path = os.path.join(current_dir, 'data')
    log.log.info('Create input data folder: %s' % input_path)
    if not os.path.exists(input_path):
        os.makedirs(input_path)
    else:
        log.log.warning('Data folder (%s) already exists' % input_path)

    delly_path = os.path.join(input_path, 'delly')
    log.log.info('Set delly data folder: %s' % delly_path)

    snowman_path = os.path.join(input_path, 'snowman')
    log.log.info('Set snowman data folder: %s' % snowman_path)

    dranger_path = os.path.join(input_path, 'dranger')
    log.log.info('Set dranger data folder: %s' % dranger_path)

    brass_path = os.path.join(input_path, 'brass')
    log.log.info('Set brass data folder: %s' % brass_path)

    mkdir_sv_data('delly', delly_path, log)
    delly_symlink = symlink_sv_data('delly', args.delly, delly_path, log)

    mkdir_sv_data('snowman', snowman_path, log)
    snowman_symlink = symlink_sv_data('snowman', args.snowman, snowman_path, log)

    mkdir_sv_data('dranger', dranger_path, log)
    dranger_symlink = symlink_sv_data('dranger', args.dranger, dranger_path, log)

    mkdir_sv_data('brass', brass_path, log)
    brass_symlink = symlink_sv_data('brass', args.brass, brass_path, log)

    merge_cmd = ['pcawg6_sv_merge_master.sh', args.run_id, dranger_symlink,
                 snowman_symlink, brass_symlink, delly_symlink]

    merge_cmd = shlex.split(' '.join(map(str, merge_cmd)))
    log.log.info('Prepare the pcawg6 merge sv command line:')
    log.log.info(' '.join(map(str, merge_cmd)))
    merge_proc = subprocess.Popen(merge_cmd)
    out = merge_proc.communicate()[0]
    code = merge_proc.returncode
    info = 'pcawg6_sv_merge_master.sh terminated, exit code: %s' % code
    if code == 0:
        log.log.info(info)
    else:
        log.log.error(info)
    concordance_path = os.path.join(output_dir, 'output',
                                    'sv_call_concordance')
    concordance_tar = os.path.join(output_dir, 'output',
                                   '%s_sv_call_concordance.tar.gz' %
                                   args.run_id)
    log.log.info('Archive results files in folder %s in file %s' %
                 (concordance_path, concordance_tar))
    tar = tarfile.open(concordance_tar, 'w:gz')
    tar.add(concordance_path, arcname='sv_call_concordance')
    tar.close()
    log.log.info('Process finished')


if __name__ == '__main__':
    main()
