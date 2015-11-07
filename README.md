Introduction
============
libgnc provides a set of API calls to encode and decode the GNC code, which stands for generation-based network coding (GNC) code. The code is one kind of erasure-correction codes (ECC). When transmitting GNC coded packets over a lossy network (e.g. when UDP is used), intermediate node can algebrically combine packets. The approach is called **network coding** [1][2]. GNC code is designed to be as a computationally efficient network code. GNC codes encode/decode against subsets of source packets called **generations**. 

The library at present supports two catagories of GNC codes: (pseudo-)random and band. Random GNC code encodes from generations that are randomly overlapped with each other whereas the band code encodes from generations that are overlapped consecutively. The "band" name comes from that its decoding matrix is banded.

Four decoders are provided in libgnc, with different performance tradeoff considerations. The generation-by-generation (GG) decoder has a linear-time complexity but exhibits higher code overhead. On the other hand, the overlap-aware (OA) decoder has optimized code overhead but exhibits higher complexity. These two decoders essentially can be used for all kinds of GNC code (i.e., not limited to the two codes currently provided in the library). The band (BD) and compact band (CBD) decoders can only be applied to the band GNC code. The decoders has optimized code overheads as well and their complexities are between GG and OA decoders. CBD decoder uses compact matrix representation and therefore has lower memory usage. BD decoder, on the other hand, employes pivoting techniques and therefore its decoding complexity is lower than that of CBD decoder. BD decoder uses full-size matrix representation for random access, as is heavily needed during pivoting. So BD decoder has higher memory usage.

For more details about GNC codes and the decoder design, please refer to [3].

Usage
============

Limitation
============

Reference
============
[1] Ahlswede, Rudolf; N. Cai, Shuo-Yen Robert Li, and Raymond Wai-Ho Yeung. "Network Information Flow". IEEE Transactions on Information Theory, IT-46 46 (4): 1204–1216, 2000.

[2] S. Li, R. Yeung, and N. Cai, "Linear Network Coding", in IEEE Transactions on Information Theory, Vol 49, No. 2, pp. 371–381, 2003.

[3] Ye Li, "Efficient Network Coding for Different Network Topologies", Queen's University PhD Thesis, Oct., 2014. Available: https://qspace.library.queensu.ca/bitstream/1974/12602/1/Li_Ye_201410_PhD.pdf

[4] Ye Li, W.-Y. Chan, and S. D. Blostein, "Network Coding with Unequal Size Overlapping Generations", in Proceedings of International Symposium on Network Coding (NetCod), pp. 161-166, Cambridge, MA, June, 2012.

[5] Ye Li, S. D. Blostein, and W.-Y. Chan, "Large File Distribution Using Efficient Generation-based Network Coding", in Proc. IEEE Globecom Workshop on Cloud Computing Systems, Networks, and Applications, Atlanta, GA, 2013.

[6] Ye Li, W.-Y. Chan, and S. D. Blostein, "Systematic Network Coding for Transmission Over Two-Hop Lossy Links", in Proc. 27th Biennial Symposium on Communications (QBSC 2014), Kingston, ON, 2014.
