Introduction
============
libgnc provides a set of API calls to encode and decode the GNC code, which stands for generation-based network coding (GNC) code. The code is one kind of erasure-correction codes (ECC). When transmitting GNC coded packets over a lossy network (e.g. when UDP is used), intermediate node can algebrically combine packets. The approach is called **network coding**, which has been a hot research topic in the past 15 years. GNC code is designed to be as a computationally efficient network code. The key idea is to encode/decode against subsets of source packets called **generations**. 

The library at present supports two GNC codes, namely (pseudo-)random and band. Random GNC code encodes from generations that are randomly overlapped with each other whereas the band code encodes from generations that are overlapped consecutively. The "band" name comes from that its decoding matrix is banded.

Three decoders are provided in libgnc, with different tradeoff considerations between performance indices. The generation-by-generation (GG) decoder has a linear-time complexity but exhibits higher code overhead. On the other hand, the overlap-aware (OA) decoder has optimized code overhead but exhibits higher complexity. These two decoders essentially can be used for all kinds of GNC code (i.e., not limited to the two codes currently provided in the library). The band decoder, however, can only be applied to the band GNC code. The decoder has optimized code overhead and its complexity is between GG and OA decoders.

For more details about GNC codes and the decoder design, please refer to [1].

Usage
============

Limitation
============

Reference
============
[1] Ye Li, "Efficient Network Coding for Different Network Topologies", Queen's University PhD Thesis, Oct., 2014. Available: https://qspace.library.queensu.ca/bitstream/1974/12602/1/Li_Ye_201410_PhD.pdf
