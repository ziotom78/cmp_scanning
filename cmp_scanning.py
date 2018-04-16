#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

from typing import List, Union
import sys
from configparser import ConfigParser, ExtendedInterpolation
from collections import namedtuple
from glob import glob
import logging as log

import click

import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt

import healpy
import numpy as np

from astropy.io import fits

DataSource = namedtuple('DataSource',
                        ['name',
                         'file_name_mask',
                         'table_hdu',
                         'time_column',
                         'theta_column',
                         'phi_column',
                         'time_factor',
                         'angle_factor',
                         'subplot'])

Configuration = namedtuple('Configuration',
                           ['nside',
                            'delta_time',
                            'number_of_frames',
                            'show_fsky',
                            'figure_width',
                            'figure_height',
                            'figure_file_name_mask',
                            'time_measure_unit',
                            'data_sources'])


def int_or_str(s: str) -> Union[int, str]:
    '''If "s" contains an integer, convert it toan integer and return it. Otherwise, return s unchanged'''

    try:
        intval = int(s)
    except ValueError:
        return s

    return intval


def read_configuration(conf_parser: ConfigParser) -> Configuration:
    ''''Build a Configuration object and initialize its fields using the parameters read by the ConfigParser object'''
    conf_sect = conf_parser['config']
    nside = conf_sect.getint('nside', 128)
    if not healpy.isnsideok(nside):
        log.error('error, wrong NSIDE (%d)', nside)
        sys.exit(1)

    delta_time = conf_sect.getfloat('delta_time')
    number_of_frames = conf_sect.getint('number_of_frames')
    show_fsky = conf_sect.getboolean('show_fsky', True)
    figure_width = conf_sect.getfloat('figure_width')
    figure_height = conf_sect.getfloat('figure_height')
    figure_file_name_mask = conf_sect.get(
        'figure_file_name_mask', 'anim%04d.png')
    time_measure_unit = conf_sect.get('time_measure_unit', '')
    data_source_names = [x.strip()
                         for x in conf_sect.get('data_sources').split(',')]

    data_sources = []  # type: List[DataSource]
    for cur_name in data_source_names:
        source_sect = conf_parser[cur_name]
        file_name_mask = source_sect.get('file_name_mask')
        table_hdu = int_or_str(source_sect.get('table_hdu'))
        time_column = int_or_str(source_sect.get('time_column'))
        theta_column = int_or_str(source_sect.get('theta_column'))
        phi_column = int_or_str(source_sect.get('phi_column'))
        time_factor = source_sect.getfloat('time_factor', 1.0)
        angle_factor = source_sect.getfloat('angle_factor', 1.0)
        subplot = source_sect.getint('subplot', 111)

        data_sources.append(DataSource(name=cur_name,
                                       file_name_mask=file_name_mask,
                                       table_hdu=table_hdu,
                                       time_column=time_column,
                                       theta_column=theta_column,
                                       phi_column=phi_column,
                                       time_factor=time_factor,
                                       angle_factor=angle_factor,
                                       subplot=subplot))

    return Configuration(nside=nside,
                         delta_time=delta_time,
                         number_of_frames=number_of_frames,
                         show_fsky=show_fsky,
                         figure_width=figure_width,
                         figure_height=figure_height,
                         figure_file_name_mask=figure_file_name_mask,
                         time_measure_unit=time_measure_unit,
                         data_sources=data_sources)


class Pointings:
    def __init__(self, time, pixidx):
        self.time = time
        self.pixidx = pixidx


def read_pointings(data_source: DataSource, file_name: str, nside: int):
    log.info('reading file "{0}"'.format(file_name))
    with fits.open(file_name) as f:
        hdu = f[data_source.table_hdu]
        time = hdu.data.field(data_source.time_column) * \
            data_source.time_factor
        theta, phi = [hdu.data.field(x) * data_source.angle_factor
                      for x in (data_source.theta_column, data_source.phi_column)]

        pixidx = healpy.ang2pix(nside, theta, phi)

    return Pointings(time=time, pixidx=pixidx)


class TodCollection:
    def __init__(self, data_source: DataSource, nside: int):
        self.data_source = data_source
        self.file_names = sorted(glob(data_source.file_name_mask))
        self.pointings = None
        self.cur_idx = 0
        self.first_time = None
        self.nside = nside
        self.map_pixels = np.zeros(healpy.nside2npix(nside))

    def get_pixidx(self, start_time: float, end_time: float):
        if (self.pointings is None) or (self.pointings.time[-1] < end_time):
            new = read_pointings(
                self.data_source, self.file_names[self.cur_idx], self.nside)

            if self.first_time is None:
                self.first_time = new.time[0]
            new.time -= self.first_time

            if self.pointings is not None:
                mask = self.pointings.time >= start_time
                self.pointings.time = np.concatenate(
                    (self.pointings.time[mask], new.time))
                self.pointings.pixidx = np.concatenate(
                    (self.pointings.pixidx[mask], new.pixidx))
            else:
                self.pointings = new

            self.cur_idx += 1

        mask = (self.pointings.time >= start_time) & (
            self.pointings.time < end_time)
        return self.pointings.pixidx[mask]


def fsky(pixels):
    '''Return the percentage of sky seen (nonzero pixels are considered to have been «seen»).'''

    ones = len(pixels[pixels > 0])
    return (100.0 * ones) / len(pixels)


@click.command()
@click.option('--save-fsky', type=str, default=None,
              help='Save the values of f_sky for each frame in a text file')
@click.argument('parameter_file')
def main(save_fsky, parameter_file):
    log.basicConfig(
        level=log.INFO, format='[%(asctime)s %(levelname)s] %(message)s')

    conf_parser = ConfigParser(interpolation=ExtendedInterpolation())
    with open(parameter_file, 'rt') as f:
        conf_parser.read_file(f)
    config = read_configuration(conf_parser)

    plt.figure(figsize=(config.figure_width, config.figure_height))
    start_time = 0

    collections = [TodCollection(data_source=x, nside=config.nside)
                   for x in config.data_sources]
    if save_fsky:
        fsky_file = open(save_fsky, 'wt')
    else:
        fsky_file = None

    for cur_frame in range(config.number_of_frames):
        plt.clf()

        if fsky_file:
            fsky_line = '{0:6d}'.format(cur_frame)

        for coll in collections:
            pixidx = coll.get_pixidx(
                start_time, start_time + config.delta_time)
            coll.map_pixels[pixidx] = 2

            cur_fsky = fsky(coll.map_pixels)
            if config.show_fsky:
                title_template = '{name} ({time:.1f}, fsky={fsky:.2f})'
            else:
                title_template = '{name} ({time:.1f})'
            map_title = title_template.format(name=coll.data_source.name,
                                              time=start_time + config.delta_time * 0.5,
                                              unit=config.time_measure_unit,
                                              fsky=cur_fsky)
            healpy.mollview(coll.map_pixels, title=map_title,
                            sub=coll.data_source.subplot, cbar=False)

            if fsky_file:
                fsky_line += '\t{0:.3f}'.format(cur_fsky)

        if fsky_file:
            fsky_file.write(fsky_line)
            fsky_file.write('\n')

        file_name = config.figure_file_name_mask % cur_frame
        plt.savefig(file_name)
        log.info('File "{0}" saved'.format(file_name))

        for coll in collections:
            coll.map_pixels[coll.map_pixels > 0] = 1

        start_time += config.delta_time

    if fsky_file:
        fsky_file.close()


if __name__ == '__main__':
    main()
