#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Fft Accum Avg P1 Rts 2048 Ettus
# Generated: Thu Nov 29 15:12:22 2018
##################################################

from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import fft
from gnuradio import gr
from gnuradio import uhd
from gnuradio.eng_option import eng_option
from gnuradio.fft import window
from gnuradio.filter import firdes
from optparse import OptionParser
import time


class fft_accum_avg_p1_rts_2048_ettus(gr.top_block):

    def __init__(self):
        gr.top_block.__init__(self, "Fft Accum Avg P1 Rts 2048 Ettus")

        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = 10e6
        self.f_stop = f_stop = 150.0e6
        self.f_start = f_start = 90.0e6
        self.vec_fft_len = vec_fft_len = 2048
        self.str_file_out = str_file_out = '/Users/josaitis/gnuradio_tutorials/fft_accum_temp/fft_accum_temp_TEST.txt'
        self.n_head = n_head = 3
        self.n_FFTs = n_FFTs = (f_stop-f_start)/samp_rate
        self.f_c = f_c = 106.5e6
        self.N_to_avg = N_to_avg = 10

        ##################################################
        # Blocks
        ##################################################
        self.uhd_usrp_source_0 = uhd.usrp_source(
        	",".join(("", "")),
        	uhd.stream_args(
        		cpu_format="fc32",
        		channels=range(1),
        	),
        )
        self.uhd_usrp_source_0.set_samp_rate(samp_rate)
        self.uhd_usrp_source_0.set_center_freq(f_c, 0)
        self.uhd_usrp_source_0.set_gain(20, 0)
        self.uhd_usrp_source_0.set_antenna('TX/RX', 0)
        self.uhd_usrp_source_0.set_auto_dc_offset(True, 0)
        self.source_to_fft_vector = blocks.stream_to_vector(gr.sizeof_gr_complex*1, vec_fft_len)
        self.fft_vxx_0 = fft.fft_vcc(vec_fft_len, True, (window.blackmanharris(vec_fft_len)), True, 1)
        self.blocks_head_0 = blocks.head(gr.sizeof_float*vec_fft_len, n_head)
        self.blocks_file_sink_0 = blocks.file_sink(gr.sizeof_float*vec_fft_len, str_file_out, True)
        self.blocks_file_sink_0.set_unbuffered(False)
        self.blocks_complex_to_mag_1 = blocks.complex_to_mag(vec_fft_len)



        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_complex_to_mag_1, 0), (self.blocks_head_0, 0))
        self.connect((self.blocks_head_0, 0), (self.blocks_file_sink_0, 0))
        self.connect((self.fft_vxx_0, 0), (self.blocks_complex_to_mag_1, 0))
        self.connect((self.source_to_fft_vector, 0), (self.fft_vxx_0, 0))
        self.connect((self.uhd_usrp_source_0, 0), (self.source_to_fft_vector, 0))

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.uhd_usrp_source_0.set_samp_rate(self.samp_rate)
        self.set_n_FFTs((self.f_stop-self.f_start)/self.samp_rate)

    def get_f_stop(self):
        return self.f_stop

    def set_f_stop(self, f_stop):
        self.f_stop = f_stop
        self.set_n_FFTs((self.f_stop-self.f_start)/self.samp_rate)

    def get_f_start(self):
        return self.f_start

    def set_f_start(self, f_start):
        self.f_start = f_start
        self.set_n_FFTs((self.f_stop-self.f_start)/self.samp_rate)

    def get_vec_fft_len(self):
        return self.vec_fft_len

    def set_vec_fft_len(self, vec_fft_len):
        self.vec_fft_len = vec_fft_len

    def get_str_file_out(self):
        return self.str_file_out

    def set_str_file_out(self, str_file_out):
        self.str_file_out = str_file_out
        self.blocks_file_sink_0.open(self.str_file_out)

    def get_n_head(self):
        return self.n_head

    def set_n_head(self, n_head):
        self.n_head = n_head
        self.blocks_head_0.set_length(self.n_head)

    def get_n_FFTs(self):
        return self.n_FFTs

    def set_n_FFTs(self, n_FFTs):
        self.n_FFTs = n_FFTs

    def get_f_c(self):
        return self.f_c

    def set_f_c(self, f_c):
        self.f_c = f_c
        self.uhd_usrp_source_0.set_center_freq(self.f_c, 0)

    def get_N_to_avg(self):
        return self.N_to_avg

    def set_N_to_avg(self, N_to_avg):
        self.N_to_avg = N_to_avg


def main(top_block_cls=fft_accum_avg_p1_rts_2048_ettus, options=None):
    if gr.enable_realtime_scheduling() != gr.RT_OK:
        print "Error: failed to enable real-time scheduling."

    tb = top_block_cls()
    tb.start()
    tb.wait()


if __name__ == '__main__':
    main()
