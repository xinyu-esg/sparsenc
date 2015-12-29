from pysnc import *
import os
import sys

def print_usage():
    print("Usage: ./ProgramName c_type d_type filename")
    print("        c_type   - code type, available options: RAND BAND")
    print("        d_type   - decoder type, available options: GG OA BD CBD")
    print("        filename - filename to be encoded and decoded.")

if len(sys.argv) != 4:
    print_usage()
    sys.exit(1)

c_type = 0
if sys.argv[1] == "RAND":
    c_type = RAND_SNC
elif sys.argv[1] == "BAND":
    c_type = BAND_SNC
else:
    print_usage()
    sys.exit(1)

d_type = 0
if sys.argv[2] == "GG":
    d_type = GG_DECODER
elif sys.argv[2] == "OA":
    d_type = OA_DECODER
elif sys.argv[2] == "BD":
    d_type = BD_DECODER
elif sys.argv[2] == "CBD":
    d_type = CBD_DECODER
else:
    print_usage()
    sys.exit(1)

filename = sys.argv[3]
if not os.path.isfile(filename):
    print("Filename is invalid.")
    sys.exit(1)

filename = filename.encode('UTF-8') # byte string literal for compatibility in Python 3.x
statinfo = os.stat(filename)
filesize = statinfo.st_size
datasize = filesize

sp = snc_parameter(datasize, 0.01, 32, 48, 1280, c_type, 0, 0)
sc = snc.snc_create_enc_context(None, sp)
snc.snc_load_file_to_context(c_char_p(filename), 0, sc)  # Load file to snc_context
decoder = snc.snc_create_decoder(sp, d_type)  # Create decoder

while not snc.snc_decoder_finished(decoder):
    pkt_p = snc.snc_generate_packet(sc)
    snc.snc_process_packet(decoder, pkt_p)

copyname = filename + b'.dec.copy'
snc.snc_recover_to_file(c_char_p(copyname), snc.snc_get_enc_context(decoder))

snc.snc_free_decoder(decoder)
snc.snc_free_enc_context(sc)
