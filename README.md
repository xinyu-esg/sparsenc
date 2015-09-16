Introduction
============
libgnc provides a set of API calls to encode and decode the GNC code, which stands for generation-based network coding (GNC) code. The code is one kind of erasure-correction codes (ECC). When transmitting GNC coded packets over a lossy network (e.g. when UDP is used), intermediate node can algebrically combine packets. The approach is called **network coding** [1][2]. GNC code is designed to be as a computationally efficient network code. GNC codes encode/decode against subsets of source packets called **generations**. 

The library at present supports two GNC codes, namely (pseudo-)random and band. Random GNC code encodes from generations that are randomly overlapped with each other whereas the band code encodes from generations that are overlapped consecutively. The "band" name comes from that its decoding matrix is banded.

Four decoders are provided in libgnc, with different performance tradeoff considerations. The generation-by-generation (GG) decoder has a linear-time complexity but exhibits higher code overhead. On the other hand, the overlap-aware (OA) decoder has optimized code overhead but exhibits higher complexity. These two decoders essentially can be used for all kinds of GNC code (i.e., not limited to the two codes currently provided in the library). The band (BD) and compact band (CBD) decoders, however, can only be applied to the band GNC code. The decoders has optimized code overheads as well and their complexities are between GG and OA decoders. CBD decoder uses compact matrix representation and therefore has lower memory usage. BD decoder, on the other hand, can employ pivoting techniques thanks to the random access capability brought by the full-size matrix representation, and hence its decoding complexity is lower than that of CBD decoder.

For more details about GNC codes and the decoder design, please refer to [3].

Usage
============

Limitation
============
An important componenet missing in libgnc is re-encoding. As a network code, GNC code allows intermediate nodes to re-encode packets such that the network throughput can be maximized. Re-encoding requires a buffer to save received GNC packets, and a scheduling algorithm to determine from which received packets to re-encode a new packet. The implementation is left as a future work. 

Reference
============
[1] Ahlswede, Rudolf; N. Cai, Shuo-Yen Robert Li, and Raymond Wai-Ho Yeung. "Network Information Flow". IEEE Transactions on Information Theory, IT-46 46 (4): 1204–1216, 2000.

[2] S. Li, R. Yeung, and N. Cai, "Linear Network Coding", in IEEE Transactions on Information Theory, Vol 49, No. 2, pp. 371–381, 2003

[3] Ye Li, "Efficient Network Coding for Different Network Topologies", Queen's University PhD Thesis, Oct., 2014. Available: https://qspace.library.queensu.ca/bitstream/1974/12602/1/Li_Ye_201410_PhD.pdf
