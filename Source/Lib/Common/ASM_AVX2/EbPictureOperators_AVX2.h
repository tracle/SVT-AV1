/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#ifndef EbPictureOperators_AVX2_h
#define EbPictureOperators_AVX2_h

#include <immintrin.h>
#include "EbDefinitions.h"
#include "EbPictureOperators_SSE2.h"

#ifdef __cplusplus
extern "C" {
#endif

uint64_t spatial_full_distortion_kernel4x_n_avx2_intrin(uint8_t *input, uint32_t input_offset,
                                                        uint32_t input_stride, uint8_t *recon,
                                                        int32_t recon_offset,
                                                        uint32_t recon_stride, uint32_t area_width,
                                                        uint32_t area_height);

uint64_t spatial_full_distortion_kernel8x_n_avx2_intrin(uint8_t *input, uint32_t input_offset,
                                                        uint32_t input_stride, uint8_t *recon,
                                                        int32_t recon_offset,
                                                        uint32_t recon_stride, uint32_t area_width,
                                                        uint32_t area_height);

uint64_t spatial_full_distortion_kernel16x_n_avx2_intrin(uint8_t *input, uint32_t input_offset,
                                                         uint32_t input_stride, uint8_t *recon,
                                                         int32_t recon_offset,
                                                         uint32_t recon_stride, uint32_t area_width,
                                                         uint32_t area_height);

uint64_t spatial_full_distortion_kernel32x_n_avx2_intrin(uint8_t *input, uint32_t input_offset,
                                                         uint32_t input_stride, uint8_t *recon,
                                                         int32_t recon_offset,
                                                         uint32_t recon_stride, uint32_t area_width,
                                                         uint32_t area_height);

uint64_t spatial_full_distortion_kernel64x_n_avx2_intrin(uint8_t *input, uint32_t input_offset,
                                                         uint32_t input_stride, uint8_t *recon,
                                                         int32_t recon_offset,
                                                         uint32_t recon_stride, uint32_t area_width,
                                                         uint32_t area_height);

uint64_t spatial_full_distortion_kernel128x_n_avx2_intrin(
    uint8_t *input, uint32_t input_offset, uint32_t input_stride, uint8_t *recon,
    int32_t recon_offset, uint32_t recon_stride, uint32_t area_width, uint32_t area_height);

#ifdef __cplusplus
}
#endif

#endif // EbPictureOperators_AVX2_h
