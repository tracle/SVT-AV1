/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/

#include "EbDefinitions.h"
#include "immintrin.h"
#include "EbUtility.h"

EB_EXTERN EB_ALIGN(16) const uint8_t filter_type[] = {
    1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4, 1, 4};

EB_EXTERN EB_ALIGN(16) const uint8_t weak_chroma_filter[2][32] = {
    {2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4,
     2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4},
    {1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
     1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2},
};

static inline void luma_weak_filter_avx2_intrin(__m256i top, __m256i curr, __m256i bottom,
                                         __m256i curr_prev, __m256i curr_next,
                                         uint8_t *ptr_denoised, uint8_t *ptr_noise) {
    __m256i top_first_half, bottom_first_half, filter_first_half, filter_second_half,
        curr_next_first_half, curr_next_second_half, weights, curr_left_mid_first_half_weight,
        curr_left_mid_first_halflo, curr_left_mid_first_halfhi, curr_prev_permutation,
        curr_permutation, curr_next_permutation, top_permutation, bottom_permutation;

    curr_prev_permutation           = _mm256_permute4x64_epi64(curr_prev, 216);
    curr_permutation                = _mm256_permute4x64_epi64(curr, 216);
    curr_left_mid_first_halflo      = _mm256_unpacklo_epi8(curr_prev_permutation, curr_permutation);
    weights                         = _mm256_loadu_si256((__m256i *)filter_type);
    curr_left_mid_first_half_weight = _mm256_maddubs_epi16(curr_left_mid_first_halflo, weights);
    curr_next_permutation           = _mm256_permute4x64_epi64(curr_next, 88);
    curr_next_first_half = _mm256_unpacklo_epi8(curr_next_permutation, _mm256_setzero_si256());
    curr_left_mid_first_halflo =
        _mm256_add_epi16(curr_next_first_half, curr_left_mid_first_half_weight);

    curr_left_mid_first_halfhi      = _mm256_unpackhi_epi8(curr_prev_permutation, curr_permutation);
    curr_left_mid_first_half_weight = _mm256_maddubs_epi16(curr_left_mid_first_halfhi, weights);
    curr_next_permutation           = _mm256_permute4x64_epi64(curr_next, 216);
    curr_next_second_half = _mm256_unpackhi_epi8(curr_next_permutation, _mm256_setzero_si256());
    curr_left_mid_first_halfhi =
        _mm256_add_epi16(curr_next_second_half, curr_left_mid_first_half_weight);

    top_permutation    = _mm256_permute4x64_epi64(top, 216);
    top_first_half     = _mm256_unpacklo_epi8(top_permutation, _mm256_setzero_si256());
    bottom_permutation = _mm256_permute4x64_epi64(bottom, 216);
    bottom_first_half  = _mm256_unpacklo_epi8(bottom_permutation, _mm256_setzero_si256());
    filter_first_half  = _mm256_adds_epi16(_mm256_adds_epi16(bottom_first_half, top_first_half),
                                          curr_left_mid_first_halflo);
    filter_first_half  = _mm256_srli_epi16(filter_first_half, 3);

    top_first_half     = _mm256_unpackhi_epi8(top_permutation, _mm256_setzero_si256());
    bottom_first_half  = _mm256_unpackhi_epi8(bottom_permutation, _mm256_setzero_si256());
    filter_second_half = _mm256_adds_epi16(_mm256_adds_epi16(bottom_first_half, top_first_half),
                                           curr_left_mid_first_halfhi);
    filter_second_half = _mm256_srli_epi16(filter_second_half, 3);

    filter_first_half =
        _mm256_permute4x64_epi64(_mm256_packus_epi16(filter_first_half, filter_second_half), 216);
    _mm256_storeu_si256((__m256i *)(ptr_denoised), filter_first_half);

    _mm256_storeu_si256((__m256i *)(ptr_noise), _mm256_subs_epu8(curr, filter_first_half));
}
static inline void luma_weak_filter_128_avx2_intrin(__m128i top, __m128i curr, __m128i bottom,
                                             __m128i curr_prev, __m128i curr_next,
                                             uint8_t *ptr_denoised, uint8_t *ptr_noise) {

    __m128i top_first_half, bottom_first_half, filter_first_half, filter_second_half,
        curr_next_first_half, curr_next_second_half, weights, curr_left_mid_first_half_weight,
        curr_left_mid_first_halflo, curr_left_mid_first_halfhi;

    curr_left_mid_first_halflo      = _mm_unpacklo_epi8(curr_prev, curr);
    weights                         = _mm_loadu_si128((__m128i *)filter_type);
    curr_left_mid_first_half_weight = _mm_maddubs_epi16(curr_left_mid_first_halflo, weights);
    curr_next_first_half            = _mm_unpacklo_epi8(curr_next, _mm_setzero_si128());
    curr_left_mid_first_halflo =
        _mm_add_epi16(curr_next_first_half, curr_left_mid_first_half_weight);

    curr_left_mid_first_halfhi      = _mm_unpackhi_epi8(curr_prev, curr);
    curr_left_mid_first_half_weight = _mm_maddubs_epi16(curr_left_mid_first_halfhi, weights);
    curr_next_second_half           = _mm_unpackhi_epi8(curr_next, _mm_setzero_si128());
    curr_left_mid_first_halfhi =
        _mm_add_epi16(curr_next_second_half, curr_left_mid_first_half_weight);

    top_first_half    = _mm_unpacklo_epi8(top, _mm_setzero_si128());
    bottom_first_half = _mm_unpacklo_epi8(bottom, _mm_setzero_si128());
    filter_first_half = _mm_adds_epi16(_mm_adds_epi16(bottom_first_half, top_first_half),
                                       curr_left_mid_first_halflo);
    filter_first_half = _mm_srli_epi16(filter_first_half, 3);

    top_first_half     = _mm_unpackhi_epi8(top, _mm_setzero_si128());
    bottom_first_half  = _mm_unpackhi_epi8(bottom, _mm_setzero_si128());
    filter_second_half = _mm_adds_epi16(_mm_adds_epi16(bottom_first_half, top_first_half),
                                        curr_left_mid_first_halfhi);
    filter_second_half = _mm_srli_epi16(filter_second_half, 3);

    filter_first_half = _mm_packus_epi16(filter_first_half, filter_second_half);
    _mm_storel_epi64((__m128i *)(ptr_denoised), filter_first_half);

    _mm_storel_epi64((__m128i *)(ptr_noise), _mm_subs_epu8(curr, filter_first_half));
}
static inline void chroma_weak_luma_strong_filter_avx2_intrin(__m256i top, __m256i curr, __m256i bottom,
                                                       __m256i curr_prev, __m256i curr_next,
                                                       __m256i top_prev, __m256i top_next,
                                                       __m256i bottom_prev, __m256i bottom_next,
                                                       uint8_t *ptr_denoised) {
    __m256i filter_first_half, filter_second_half, curr_next_first_half, curr_next_second_half,
        weights, curr_left_mid_first_half_weight, curr_left_mid_first_halflo,
        curr_left_mid_first_halfhi, curr_prev_permutation, curr_permutation, curr_next_permutation,
        top_permutation, bottom_permutation, top_prev_permutation, top_left_mid_first_halflo,
        top_left_mid_first_half_weight, top_next_first_half, top_next_permutation,
        top_left_mid_first_halfhi, top_next_second_half, bottom_prev_permutation,
        bottom_left_mid_first_halflo, bottom_left_mid_first_half_weight, bottom_next_permutation,
        bottom_next_first_half, bottom_left_mid_first_halfhi, bottom_next_second_half;

    //  Curr
    curr_prev_permutation           = _mm256_permute4x64_epi64(curr_prev, 216);
    curr_permutation                = _mm256_permute4x64_epi64(curr, 216);
    curr_left_mid_first_halflo      = _mm256_unpacklo_epi8(curr_prev_permutation, curr_permutation);
    weights                         = _mm256_loadu_si256((__m256i *)weak_chroma_filter[0]);
    curr_left_mid_first_half_weight = _mm256_maddubs_epi16(curr_left_mid_first_halflo, weights);
    curr_next_permutation           = _mm256_permute4x64_epi64(curr_next, 88);
    curr_next_first_half = _mm256_unpacklo_epi8(curr_next_permutation, _mm256_setzero_si256());
    curr_next_first_half = _mm256_slli_epi16(curr_next_first_half, 1);
    curr_left_mid_first_halflo =
        _mm256_add_epi16(curr_next_first_half, curr_left_mid_first_half_weight);

    curr_left_mid_first_halfhi      = _mm256_unpackhi_epi8(curr_prev_permutation, curr_permutation);
    curr_left_mid_first_half_weight = _mm256_maddubs_epi16(curr_left_mid_first_halfhi, weights);
    curr_next_permutation           = _mm256_permute4x64_epi64(curr_next, 216);
    curr_next_second_half = _mm256_unpackhi_epi8(curr_next_permutation, _mm256_setzero_si256());
    curr_next_second_half = _mm256_slli_epi16(curr_next_second_half, 1);
    curr_left_mid_first_halfhi =
        _mm256_add_epi16(curr_next_second_half, curr_left_mid_first_half_weight);

    // Top
    top_prev_permutation           = _mm256_permute4x64_epi64(top_prev, 216);
    top_permutation                = _mm256_permute4x64_epi64(top, 216);
    top_left_mid_first_halflo      = _mm256_unpacklo_epi8(top_prev_permutation, top_permutation);
    weights                        = _mm256_loadu_si256((__m256i *)weak_chroma_filter[1]);
    top_left_mid_first_half_weight = _mm256_maddubs_epi16(top_left_mid_first_halflo, weights);
    top_next_permutation           = _mm256_permute4x64_epi64(top_next, 88);
    top_next_first_half = _mm256_unpacklo_epi8(top_next_permutation, _mm256_setzero_si256());
    top_left_mid_first_halflo =
        _mm256_add_epi16(top_next_first_half, top_left_mid_first_half_weight);

    top_left_mid_first_halfhi      = _mm256_unpackhi_epi8(top_prev_permutation, top_permutation);
    top_left_mid_first_half_weight = _mm256_maddubs_epi16(top_left_mid_first_halfhi, weights);
    top_next_permutation           = _mm256_permute4x64_epi64(top_next, 216);
    top_next_second_half = _mm256_unpackhi_epi8(top_next_permutation, _mm256_setzero_si256());
    top_left_mid_first_halfhi =
        _mm256_add_epi16(top_next_second_half, top_left_mid_first_half_weight);

    // Bottom
    bottom_prev_permutation = _mm256_permute4x64_epi64(bottom_prev, 216);
    bottom_permutation      = _mm256_permute4x64_epi64(bottom, 216);
    bottom_left_mid_first_halflo =
        _mm256_unpacklo_epi8(bottom_prev_permutation, bottom_permutation);
    weights                           = _mm256_loadu_si256((__m256i *)weak_chroma_filter[1]);
    bottom_left_mid_first_half_weight = _mm256_maddubs_epi16(bottom_left_mid_first_halflo, weights);
    bottom_next_permutation           = _mm256_permute4x64_epi64(bottom_next, 88);
    bottom_next_first_half = _mm256_unpacklo_epi8(bottom_next_permutation, _mm256_setzero_si256());
    bottom_left_mid_first_halflo =
        _mm256_add_epi16(bottom_next_first_half, bottom_left_mid_first_half_weight);

    bottom_left_mid_first_halfhi =
        _mm256_unpackhi_epi8(bottom_prev_permutation, bottom_permutation);
    bottom_left_mid_first_half_weight = _mm256_maddubs_epi16(bottom_left_mid_first_halfhi, weights);
    bottom_next_permutation           = _mm256_permute4x64_epi64(bottom_next, 216);
    bottom_next_second_half = _mm256_unpackhi_epi8(bottom_next_permutation, _mm256_setzero_si256());
    bottom_left_mid_first_halfhi =
        _mm256_add_epi16(bottom_next_second_half, bottom_left_mid_first_half_weight);

    filter_first_half = _mm256_adds_epi16(
        _mm256_adds_epi16(bottom_left_mid_first_halflo, top_left_mid_first_halflo),
        curr_left_mid_first_halflo);
    filter_first_half  = _mm256_srli_epi16(filter_first_half, 4);
    filter_second_half = _mm256_adds_epi16(
        _mm256_adds_epi16(bottom_left_mid_first_halfhi, top_left_mid_first_halfhi),
        curr_left_mid_first_halfhi);
    filter_second_half = _mm256_srli_epi16(filter_second_half, 4);

    filter_first_half =
        _mm256_permute4x64_epi64(_mm256_packus_epi16(filter_first_half, filter_second_half), 216);
    _mm256_storeu_si256((__m256i *)(ptr_denoised), filter_first_half);
}

static inline void chroma_weak_luma_strong_filter_128_avx2_intrin(__m128i top, __m128i curr,
                                                           __m128i bottom, __m128i curr_prev,
                                                           __m128i curr_next, __m128i top_prev,
                                                           __m128i top_next, __m128i bottom_prev,
                                                           __m128i  bottom_next,
                                                           uint8_t *ptr_denoised) {
    __m128i filter_first_half, filter_second_half, curr_next_first_half, curr_next_second_half,
        weights, curr_left_mid_first_half_weight, curr_left_mid_first_halflo,
        curr_left_mid_first_halfhi, top_left_mid_first_halflo, top_left_mid_first_half_weight,
        top_next_first_half, top_left_mid_first_halfhi, top_next_second_half,
        bottom_left_mid_first_halflo, bottom_left_mid_first_half_weight, bottom_next_first_half,
        bottom_left_mid_first_halfhi, bottom_next_second_half;

    //  Curr
    curr_left_mid_first_halflo      = _mm_unpacklo_epi8(curr_prev, curr);
    weights                         = _mm_loadu_si128((__m128i *)weak_chroma_filter[0]);
    curr_left_mid_first_half_weight = _mm_maddubs_epi16(curr_left_mid_first_halflo, weights);
    curr_next_first_half            = _mm_unpacklo_epi8(curr_next, _mm_setzero_si128());
    curr_next_first_half            = _mm_slli_epi16(curr_next_first_half, 1);
    curr_left_mid_first_halflo =
        _mm_add_epi16(curr_next_first_half, curr_left_mid_first_half_weight);

    curr_left_mid_first_halfhi      = _mm_unpackhi_epi8(curr_prev, curr);
    curr_left_mid_first_half_weight = _mm_maddubs_epi16(curr_left_mid_first_halfhi, weights);
    curr_next_second_half           = _mm_unpackhi_epi8(curr_next, _mm_setzero_si128());
    curr_next_second_half           = _mm_slli_epi16(curr_next_second_half, 1);
    curr_left_mid_first_halfhi =
        _mm_add_epi16(curr_next_second_half, curr_left_mid_first_half_weight);

    // Top
    top_left_mid_first_halflo      = _mm_unpacklo_epi8(top_prev, top);
    weights                        = _mm_loadu_si128((__m128i *)weak_chroma_filter[1]);
    top_left_mid_first_half_weight = _mm_maddubs_epi16(top_left_mid_first_halflo, weights);
    top_next_first_half            = _mm_unpacklo_epi8(top_next, _mm_setzero_si128());
    top_left_mid_first_halflo = _mm_add_epi16(top_next_first_half, top_left_mid_first_half_weight);

    top_left_mid_first_halfhi      = _mm_unpackhi_epi8(top_prev, top);
    top_left_mid_first_half_weight = _mm_maddubs_epi16(top_left_mid_first_halfhi, weights);
    top_next_second_half           = _mm_unpackhi_epi8(top_next, _mm_setzero_si128());
    top_left_mid_first_halfhi = _mm_add_epi16(top_next_second_half, top_left_mid_first_half_weight);

    // Bottom
    bottom_left_mid_first_halflo      = _mm_unpacklo_epi8(bottom_prev, bottom);
    weights                           = _mm_loadu_si128((__m128i *)weak_chroma_filter[1]);
    bottom_left_mid_first_half_weight = _mm_maddubs_epi16(bottom_left_mid_first_halflo, weights);
    bottom_next_first_half            = _mm_unpacklo_epi8(bottom_next, _mm_setzero_si128());
    bottom_left_mid_first_halflo =
        _mm_add_epi16(bottom_next_first_half, bottom_left_mid_first_half_weight);

    bottom_left_mid_first_halfhi      = _mm_unpackhi_epi8(bottom_prev, bottom);
    bottom_left_mid_first_half_weight = _mm_maddubs_epi16(bottom_left_mid_first_halfhi, weights);
    bottom_next_second_half           = _mm_unpackhi_epi8(bottom_next, _mm_setzero_si128());
    bottom_left_mid_first_halfhi =
        _mm_add_epi16(bottom_next_second_half, bottom_left_mid_first_half_weight);

    filter_first_half =
        _mm_adds_epi16(_mm_adds_epi16(bottom_left_mid_first_halflo, top_left_mid_first_halflo),
                       curr_left_mid_first_halflo);
    filter_first_half = _mm_srli_epi16(filter_first_half, 4);
    filter_second_half =
        _mm_adds_epi16(_mm_adds_epi16(bottom_left_mid_first_halfhi, top_left_mid_first_halfhi),
                       curr_left_mid_first_halfhi);
    filter_second_half = _mm_srli_epi16(filter_second_half, 4);

    filter_first_half = _mm_packus_epi16(filter_first_half, filter_second_half);
    _mm_storel_epi64((__m128i *)(ptr_denoised), filter_first_half);
}

static inline void chroma_strong_avx2_intrin(__m256i top, __m256i curr, __m256i bottom, __m256i curr_prev,
                                      __m256i curr_next, __m256i top_prev, __m256i top_next,
                                      __m256i bottom_prev, __m256i bottom_next,
                                      uint8_t *ptr_denoised) {
    __m256i curr_left_mid_first_halflo, curr_left_mid_first_halfhi, curr_prev_permutation,
        curr_permutation, curr_next_permutation, top_permutation, top_prev_permutation,
        top_left_mid_first_halflo, top_next_permutation, top_left_mid_first_halfhi,
        bottom_permutation, bottom_prev_permutation, bottom_left_mid_first_halflo,
        bottom_next_permutation, bottom_left_mid_first_halfhi;

    curr_prev_permutation = _mm256_permute4x64_epi64(curr_prev, 216);
    curr_permutation      = _mm256_permute4x64_epi64(curr, 216);
    curr_next_permutation = _mm256_permute4x64_epi64(curr_next, 216);

    curr_left_mid_first_halflo =
        _mm256_add_epi16(_mm256_unpacklo_epi8(curr_permutation, _mm256_setzero_si256()),
                         _mm256_unpacklo_epi8(curr_prev_permutation, _mm256_setzero_si256()));
    curr_left_mid_first_halflo =
        _mm256_add_epi16(_mm256_unpacklo_epi8(curr_next_permutation, _mm256_setzero_si256()),
                         curr_left_mid_first_halflo);

    curr_left_mid_first_halfhi =
        _mm256_add_epi16(_mm256_unpackhi_epi8(curr_permutation, _mm256_setzero_si256()),
                         _mm256_unpackhi_epi8(curr_prev_permutation, _mm256_setzero_si256()));
    curr_left_mid_first_halfhi =
        _mm256_add_epi16(_mm256_unpackhi_epi8(curr_next_permutation, _mm256_setzero_si256()),
                         curr_left_mid_first_halfhi);

    top_prev_permutation = _mm256_permute4x64_epi64(top_prev, 216);
    top_permutation      = _mm256_permute4x64_epi64(top, 216);
    top_next_permutation = _mm256_permute4x64_epi64(top_next, 216);

    top_left_mid_first_halflo =
        _mm256_add_epi16(_mm256_unpacklo_epi8(top_permutation, _mm256_setzero_si256()),
                         _mm256_unpacklo_epi8(top_prev_permutation, _mm256_setzero_si256()));
    top_left_mid_first_halflo =
        _mm256_add_epi16(_mm256_unpacklo_epi8(top_next_permutation, _mm256_setzero_si256()),
                         top_left_mid_first_halflo);

    top_left_mid_first_halfhi =
        _mm256_add_epi16(_mm256_unpackhi_epi8(top_permutation, _mm256_setzero_si256()),
                         _mm256_unpackhi_epi8(top_prev_permutation, _mm256_setzero_si256()));
    top_left_mid_first_halfhi =
        _mm256_add_epi16(_mm256_unpackhi_epi8(top_next_permutation, _mm256_setzero_si256()),
                         top_left_mid_first_halfhi);

    bottom_prev_permutation = _mm256_permute4x64_epi64(bottom_prev, 216);
    bottom_permutation      = _mm256_permute4x64_epi64(bottom, 216);
    bottom_next_permutation = _mm256_permute4x64_epi64(bottom_next, 216);

    bottom_left_mid_first_halflo =
        _mm256_add_epi16(_mm256_unpacklo_epi8(bottom_permutation, _mm256_setzero_si256()),
                         _mm256_unpacklo_epi8(bottom_prev_permutation, _mm256_setzero_si256()));
    bottom_left_mid_first_halflo =
        _mm256_add_epi16(_mm256_unpacklo_epi8(bottom_next_permutation, _mm256_setzero_si256()),
                         bottom_left_mid_first_halflo);

    bottom_left_mid_first_halfhi =
        _mm256_add_epi16(_mm256_unpackhi_epi8(bottom_permutation, _mm256_setzero_si256()),
                         _mm256_unpackhi_epi8(bottom_prev_permutation, _mm256_setzero_si256()));
    bottom_left_mid_first_halfhi =
        _mm256_add_epi16(_mm256_unpackhi_epi8(bottom_next_permutation, _mm256_setzero_si256()),
                         bottom_left_mid_first_halfhi);

    curr_left_mid_first_halflo =
        _mm256_add_epi16(_mm256_add_epi16(curr_left_mid_first_halflo, top_left_mid_first_halflo),
                         bottom_left_mid_first_halflo);
    curr_left_mid_first_halfhi =
        _mm256_add_epi16(_mm256_add_epi16(curr_left_mid_first_halfhi, top_left_mid_first_halfhi),
                         bottom_left_mid_first_halfhi);

    top_left_mid_first_halflo =
        _mm256_unpacklo_epi16(curr_left_mid_first_halflo, _mm256_setzero_si256());
    top_left_mid_first_halflo =
        _mm256_mullo_epi32(top_left_mid_first_halflo, _mm256_set1_epi32(7282));
    top_left_mid_first_halflo = _mm256_srli_epi32(top_left_mid_first_halflo, 16);
    bottom_left_mid_first_halflo =
        _mm256_unpackhi_epi16(curr_left_mid_first_halflo, _mm256_setzero_si256());
    bottom_left_mid_first_halflo =
        _mm256_mullo_epi32(bottom_left_mid_first_halflo, _mm256_set1_epi32(7282));
    bottom_left_mid_first_halflo = _mm256_srli_epi32(bottom_left_mid_first_halflo, 16);
    curr_left_mid_first_halflo =
        _mm256_packus_epi32(top_left_mid_first_halflo, bottom_left_mid_first_halflo);

    curr_left_mid_first_halflo = _mm256_insertf128_si256(
        _mm256_setzero_si256(),
        _mm_packus_epi16(_mm256_extracti128_si256(curr_left_mid_first_halflo, 0),
                         _mm256_extracti128_si256(curr_left_mid_first_halflo, 1)),
        0);

    top_left_mid_first_halfhi =
        _mm256_unpacklo_epi16(curr_left_mid_first_halfhi, _mm256_setzero_si256());
    top_left_mid_first_halfhi =
        _mm256_mullo_epi32(top_left_mid_first_halfhi, _mm256_set1_epi32(7282));
    top_left_mid_first_halfhi = _mm256_srli_epi32(top_left_mid_first_halfhi, 16);

    bottom_left_mid_first_halfhi =
        _mm256_unpackhi_epi16(curr_left_mid_first_halfhi, _mm256_setzero_si256());
    bottom_left_mid_first_halfhi =
        _mm256_mullo_epi32(bottom_left_mid_first_halfhi, _mm256_set1_epi32(7282));
    bottom_left_mid_first_halfhi = _mm256_srli_epi32(bottom_left_mid_first_halfhi, 16);
    curr_left_mid_first_halfhi =
        _mm256_packus_epi32(top_left_mid_first_halfhi, bottom_left_mid_first_halfhi);

    curr_left_mid_first_halflo = _mm256_insertf128_si256(
        curr_left_mid_first_halflo,
        _mm_packus_epi16(_mm256_extracti128_si256(curr_left_mid_first_halfhi, 0),
                         _mm256_extracti128_si256(curr_left_mid_first_halfhi, 1)),
        1);
    _mm256_storeu_si256((__m256i *)(ptr_denoised), curr_left_mid_first_halflo);
}

static inline void chroma_strong_128_avx2_intrin(__m128i top, __m128i curr, __m128i bottom,
                                          __m128i curr_prev, __m128i curr_next, __m128i top_prev,
                                          __m128i top_next, __m128i bottom_prev,
                                          __m128i bottom_next, uint8_t *ptr_denoised) {
    __m128i curr_left_mid_first_halflo, curr_left_mid_first_halfhi, top_left_mid_first_halflo,
        top_left_mid_first_halfhi, bottom_left_mid_first_halflo, bottom_left_mid_first_halfhi;

    curr_left_mid_first_halflo = _mm_add_epi16(_mm_unpacklo_epi8(curr, _mm_setzero_si128()),
                                               _mm_unpacklo_epi8(curr_prev, _mm_setzero_si128()));
    curr_left_mid_first_halflo = _mm_add_epi16(_mm_unpacklo_epi8(curr_next, _mm_setzero_si128()),
                                               curr_left_mid_first_halflo);

    curr_left_mid_first_halfhi = _mm_add_epi16(_mm_unpackhi_epi8(curr, _mm_setzero_si128()),
                                               _mm_unpackhi_epi8(curr_prev, _mm_setzero_si128()));
    curr_left_mid_first_halfhi = _mm_add_epi16(_mm_unpackhi_epi8(curr_next, _mm_setzero_si128()),
                                               curr_left_mid_first_halfhi);

    top_left_mid_first_halflo = _mm_add_epi16(_mm_unpacklo_epi8(top, _mm_setzero_si128()),
                                              _mm_unpacklo_epi8(top_prev, _mm_setzero_si128()));
    top_left_mid_first_halflo =
        _mm_add_epi16(_mm_unpacklo_epi8(top_next, _mm_setzero_si128()), top_left_mid_first_halflo);

    top_left_mid_first_halfhi = _mm_add_epi16(_mm_unpackhi_epi8(top, _mm_setzero_si128()),
                                              _mm_unpackhi_epi8(top_prev, _mm_setzero_si128()));
    top_left_mid_first_halfhi =
        _mm_add_epi16(_mm_unpackhi_epi8(top_next, _mm_setzero_si128()), top_left_mid_first_halfhi);

    bottom_left_mid_first_halflo =
        _mm_add_epi16(_mm_unpacklo_epi8(bottom, _mm_setzero_si128()),
                      _mm_unpacklo_epi8(bottom_prev, _mm_setzero_si128()));
    bottom_left_mid_first_halflo = _mm_add_epi16(
        _mm_unpacklo_epi8(bottom_next, _mm_setzero_si128()), bottom_left_mid_first_halflo);

    bottom_left_mid_first_halfhi =
        _mm_add_epi16(_mm_unpackhi_epi8(bottom, _mm_setzero_si128()),
                      _mm_unpackhi_epi8(bottom_prev, _mm_setzero_si128()));
    bottom_left_mid_first_halfhi = _mm_add_epi16(
        _mm_unpackhi_epi8(bottom_next, _mm_setzero_si128()), bottom_left_mid_first_halfhi);

    curr_left_mid_first_halflo =
        _mm_add_epi16(_mm_add_epi16(curr_left_mid_first_halflo, top_left_mid_first_halflo),
                      bottom_left_mid_first_halflo);
    curr_left_mid_first_halfhi =
        _mm_add_epi16(_mm_add_epi16(curr_left_mid_first_halfhi, top_left_mid_first_halfhi),
                      bottom_left_mid_first_halfhi);

    top_left_mid_first_halflo = _mm_unpacklo_epi16(curr_left_mid_first_halflo, _mm_setzero_si128());
    top_left_mid_first_halflo = _mm_mullo_epi32(top_left_mid_first_halflo, _mm_set1_epi32(7282));
    top_left_mid_first_halflo = _mm_srli_epi32(top_left_mid_first_halflo, 16);
    bottom_left_mid_first_halflo =
        _mm_unpackhi_epi16(curr_left_mid_first_halflo, _mm_setzero_si128());
    bottom_left_mid_first_halflo =
        _mm_mullo_epi32(bottom_left_mid_first_halflo, _mm_set1_epi32(7282));
    bottom_left_mid_first_halflo = _mm_srli_epi32(bottom_left_mid_first_halflo, 16);
    curr_left_mid_first_halflo =
        _mm_packus_epi32(top_left_mid_first_halflo, bottom_left_mid_first_halflo);

    curr_left_mid_first_halflo =
        _mm_packus_epi16(curr_left_mid_first_halflo, curr_left_mid_first_halflo);

    _mm_storel_epi64((__m128i *)(ptr_denoised), curr_left_mid_first_halflo);
}
