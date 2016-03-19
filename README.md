Introduction
============
This project provides a set of APIs to encode and decode **sparse** **n**etwork **c**odes. Network coding [1][2] is known to be able to improve network throughput and combat packet loss by algebrically combining packets at intermediate nodes. This library provides several computationally efficient network codes to overcome the high complexity issue of the random linear network coding (RLNC). The basic idea is to encode against subsets of source packets. In early development of this library, the subsets are referred to as generations [3][4][5]. In current version, the subsets are referred to as **subgenerations** because in network coding literature the term "generation" is sometimes used to refer to the whole set of source packets to be coded by RLNC. Using subgenerations is to 1) avoid confusion and 2) to emphasize the subset concept. The RLNC is a special case that the number of subgenerations is one.

The library at present supports two major catagories of sparse network codes (SNC): random and band codes. Random SNC encodes from subgenerations that are randomly overlapped whereas the band code encodes from subgenerations that are overlapped consecutively (in terms of packet indices). The _band_ name comes from that its decoding matrix is a band matrix.

Four decoders with different performance tradeoff considerations are provided in the library. The (sub)generation-by-(sub)generation (GG) decoder has a linear-time complexity but exhibits higher code overhead. On the other hand, the overlap-aware (OA) decoder has optimized code overhead but exhibits higher complexity. These two decoders essentially can be used for all kinds of subgeneration-based code (i.e., not limited to the two codes currently provided in the library). The band (BD) and compact band (CBD) decoders can only be applied to the band code. The decoders have optimized code overheads as well and their complexities are between GG and OA decoders. CBD decoder uses compact matrix representation and therefore has lower memory usage. BD decoder, on the other hand, employes pivoting techniques and therefore its decoding complexity is lower. BD decoder uses full-size matrix representation for random access, as is heavily needed during pivoting. So BD decoder has higher memory usage.

The construction of original RLNC is can be done by constructing a band code with one subgeneration. The decoding of RLNC can be done using CBD decoder.

Systematic coding is also supported, which generates coded packets from each subgeneration only after each source packet therein is sent once. The decoding cost can be significantly reduced when the code is used in networks with low packet loss rate. The price to pay is higher overhead if the number of subgenerations is greater than 1. It is recommended to use systematic coding for RLNC in many scenarios (e.g., [6]).

For more details about subgeneration-based codes and the decoder design, please refer to [3].

Usage
============
The library is available as a shared library which is compiled by

```shell
$ make libsparsenc.so
```

Accessing API is via `include/sparsenc.h`. 

To test codes and decoders, run

```shell
$ make sncDecoders
```

and test using

```shell
usage: ./sncDecoder code_t dec_t datasize pcrate size_b size_g size_p bpc bnc sys
                       code_t   - RAND, BAND, WINDWRAP
                       dec_t    - GG, OA, BD, CBD
                       datasize - Number of bytes
                       pcrate   - Precode rate (percentage of check packets)
                       size_b   - Subgeneration distance
                       size_g   - Subgeneration size
                       size_p   - Packet size in bytes
                       bpc      - Use binary precode (0 or 1)
                       bnc      - Use binary network code (0 or 1)
                       sys      - Systematic code (0 or 1)
```

Please see `makefile` for other available examples.

Limitation
============
The library only supports coding against a given block of source packets, i.e., a *generation* of packets as termed in the network coding literature. Sliding-window mode is not supported.

Reference
============
[1] Ahlswede, Rudolf; N. Cai, Shuo-Yen Robert Li, and Raymond Wai-Ho Yeung. "Network Information Flow". IEEE Transactions on Information Theory, IT-46 46 (4): 1204–1216, 2000.

[2] S. Li, R. Yeung, and N. Cai, "Linear Network Coding", in IEEE Transactions on Information Theory, Vol 49, No. 2, pp. 371–381, 2003.

[3] Ye Li, "Efficient Network Coding for Different Network Topologies", Queen's University PhD Thesis, Oct., 2014. Available: https://qspace.library.queensu.ca/bitstream/1974/12602/1/Li_Ye_201410_PhD.pdf

[4] Ye Li, W.-Y. Chan, and S. D. Blostein, "Network Coding with Unequal Size Overlapping Generations", in Proceedings of International Symposium on Network Coding (NetCod), pp. 161-166, Cambridge, MA, June, 2012.

[5] Ye Li, S. D. Blostein, and W.-Y. Chan, "Large File Distribution Using Efficient Generation-based Network Coding", in Proc. IEEE Globecom Workshop on Cloud Computing Systems, Networks, and Applications, Atlanta, GA, 2013.

[6] Ye Li, W.-Y. Chan, and S. D. Blostein, "Systematic Network Coding for Transmission Over Two-Hop Lossy Links", in Proc. 27th Biennial Symposium on Communications (QBSC 2014), Kingston, ON, 2014. 
