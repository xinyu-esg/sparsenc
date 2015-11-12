Introduction
============
libslnc provides a set of API calls to encode and decode sparse linear network codes. Linear network coding [1][2] is able to improve network throughput and combat packet loss by algebrically combining packets at intermediate packets. This library provides several computationally efficient network codes to overcome the high complexity issue of random linear network coding. The basic idea is to encode/decode against subsets of source packets called **subgenerations**. 

The library at present supports two catagories of codes: (pseudo-)random and band. Random SLNC code encodes from subgenerations that are randomly overlapped whereas the band code encodes from subgenerations that are overlapped consecutively. The "band" name comes from that its decoding matrix is banded.

Four decoders are provided, with different performance tradeoff considerations. The (sub)generation-by-(sub)generation (GG) decoder has a linear-time complexity but exhibits higher code overhead. On the other hand, the overlap-aware (OA) decoder has optimized code overhead but exhibits higher complexity. These two decoders essentially can be used for all kinds of subgeneration-based code (i.e., not limited to the two codes currently provided in the library). The band (BD) and compact band (CBD) decoders can only be applied to the band code. The decoders have optimized code overheads as well and their complexities are between GG and OA decoders. CBD decoder uses compact matrix representation and therefore has lower memory usage. BD decoder, on the other hand, employes pivoting techniques and therefore its decoding complexity is lower. BD decoder uses full-size matrix representation for random access, as is heavily needed during pivoting. So BD decoder has higher memory usage.

For more details about subgeneration-based codes and the decoder design, please refer to [3].

Usage
============

Limitation
============
The library only supports coding against a given block of source packets, i.e., a *generation* of packets as termed in many network coding papers. Sliding-window mode is not supported.

Reference
============
[1] Ahlswede, Rudolf; N. Cai, Shuo-Yen Robert Li, and Raymond Wai-Ho Yeung. "Network Information Flow". IEEE Transactions on Information Theory, IT-46 46 (4): 1204–1216, 2000.

[2] S. Li, R. Yeung, and N. Cai, "Linear Network Coding", in IEEE Transactions on Information Theory, Vol 49, No. 2, pp. 371–381, 2003.

[3] Ye Li, "Efficient Network Coding for Different Network Topologies", Queen's University PhD Thesis, Oct., 2014. Available: https://qspace.library.queensu.ca/bitstream/1974/12602/1/Li_Ye_201410_PhD.pdf

[4] Ye Li, W.-Y. Chan, and S. D. Blostein, "Network Coding with Unequal Size Overlapping Generations", in Proceedings of International Symposium on Network Coding (NetCod), pp. 161-166, Cambridge, MA, June, 2012.

[5] Ye Li, S. D. Blostein, and W.-Y. Chan, "Large File Distribution Using Efficient Generation-based Network Coding", in Proc. IEEE Globecom Workshop on Cloud Computing Systems, Networks, and Applications, Atlanta, GA, 2013.

[6] Ye Li, W.-Y. Chan, and S. D. Blostein, "Systematic Network Coding for Transmission Over Two-Hop Lossy Links", in Proc. 27th Biennial Symposium on Communications (QBSC 2014), Kingston, ON, 2014. 
