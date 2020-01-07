/*
* Copyright(c) 2019 Intel Corporation
* SPDX - License - Identifier: BSD - 2 - Clause - Patent
*/
#include "EbEncHandle.h"
#include "EbPictureControlSet.h"
#include "EbSequenceControlSet.h"
#include "EbMotionEstimationResults.h"
#include "EbInitialRateControlProcess.h"
#include "EbInitialRateControlResults.h"
#include "EbMotionEstimationContext.h"
#include "EbUtility.h"
#include "EbReferenceObject.h"
#include "EbResize.h"
#include "EbTransforms.h"
#include "aom_dsp_rtcd.h"

/**************************************
 * Context
 **************************************/
typedef struct InitialRateControlContext {
    EbFifo *motion_estimation_results_input_fifo_ptr;
    EbFifo *initialrate_control_results_output_fifo_ptr;
} InitialRateControlContext;

/**************************************
* Macros
**************************************/
#define PAN_SB_PERCENTAGE 75
#define LOW_AMPLITUDE_TH 16

static void eb_get_mv(PictureParentControlSet *pcs_ptr, uint32_t sb_index, int32_t *x_current_mv,
            int32_t *y_current_mv) {
    uint32_t me_candidate_index;

    const MeSbResults *me_results       = pcs_ptr->me_results[sb_index];
    uint8_t            total_me_cnt     = me_results->total_me_candidate_index[0];
    const MeCandidate *me_block_results = me_results->me_candidate[0];
    for (me_candidate_index = 0; me_candidate_index < total_me_cnt; me_candidate_index++) {
        if (me_block_results->direction == UNI_PRED_LIST_0) {
            *x_current_mv = me_results->me_mv_array[0][0].x_mv;
            *y_current_mv = me_results->me_mv_array[0][0].y_mv;
            break;
        }
    }
}
EbBool check_mv_for_pan_high_amp(uint32_t hierarchical_levels, uint32_t temporal_layer_index,
                                 int32_t *x_current_mv, int32_t *x_candidate_mv) {
    if (*x_current_mv * *x_candidate_mv >
            0 // both negative or both positives and both different than 0 i.e. same direction and non Stationary)
        && ABS(*x_current_mv) >=
               global_motion_threshold[hierarchical_levels][temporal_layer_index] // high amplitude
        && ABS(*x_candidate_mv) >=
               global_motion_threshold[hierarchical_levels][temporal_layer_index] // high amplitude
        && ABS(*x_current_mv - *x_candidate_mv) < LOW_AMPLITUDE_TH) { // close amplitude

        return (EB_TRUE);

    }

    else
        return (EB_FALSE);
}

void GetMeDist(
    PictureParentControlSet    *picture_control_set_ptr,
    uint32_t                         sb_index,
    uint32_t                      *distortion)
{
    *distortion = (uint32_t)picture_control_set_ptr->me_results[sb_index]->me_candidate[0][0].distortion;
}

#if CUTREE_LA
MeCandidate* GetMbMv(
    PictureParentControlSet    *picture_control_set_ptr,
    uint32_t                    sb_index,
    uint32_t                    rf_idx,
    uint32_t                    me_mb_offset,
    int32_t                    *xCurrentMv,
    int32_t                    *yCurrentMv)
{
    uint32_t             meCandidateIndex;

    const MeLcuResults *me_results = picture_control_set_ptr->me_results[sb_index];
    uint8_t total_me_cnt = me_results->total_me_candidate_index[me_mb_offset];
    MeCandidate *me_block_results = me_results->me_candidate[me_mb_offset];
    for (meCandidateIndex = 0; meCandidateIndex < total_me_cnt; meCandidateIndex++) {
        if (me_block_results->direction == UNI_PRED_LIST_0) {
            *xCurrentMv = me_results->me_mv_array[me_mb_offset][rf_idx].x_mv;
            *yCurrentMv = me_results->me_mv_array[me_mb_offset][rf_idx].y_mv;
            break;
        }
    }
    return me_block_results;
}
#endif

EbBool check_mv_for_tilt_high_amp(uint32_t hierarchical_levels, uint32_t temporal_layer_index,
                                  int32_t *y_current_mv, int32_t *y_candidate_mv) {
    if (*y_current_mv * *y_candidate_mv >
            0 // both negative or both positives and both different than 0 i.e. same direction and non Stationary)
        && ABS(*y_current_mv) >=
               global_motion_threshold[hierarchical_levels][temporal_layer_index] // high amplitude
        && ABS(*y_candidate_mv) >=
               global_motion_threshold[hierarchical_levels][temporal_layer_index] // high amplitude
        && ABS(*y_current_mv - *y_candidate_mv) < LOW_AMPLITUDE_TH) { // close amplitude

        return (EB_TRUE);
    }

    else
        return (EB_FALSE);
}

EbBool check_mv_for_pan(uint32_t hierarchical_levels, uint32_t temporal_layer_index,
                        int32_t *x_current_mv, int32_t *y_current_mv, int32_t *x_candidate_mv,
                        int32_t *y_candidate_mv) {
    if (*y_current_mv < LOW_AMPLITUDE_TH &&
        *y_candidate_mv<
            LOW_AMPLITUDE_TH && * x_current_mv * *
            x_candidate_mv > 0 // both negative or both positives and both different than 0 i.e. same direction and non Stationary)
        && ABS(*x_current_mv) >=
               global_motion_threshold[hierarchical_levels][temporal_layer_index] // high amplitude
        && ABS(*x_candidate_mv) >=
               global_motion_threshold[hierarchical_levels][temporal_layer_index] // high amplitude
        && ABS(*x_current_mv - *x_candidate_mv) < LOW_AMPLITUDE_TH) { // close amplitude

        return (EB_TRUE);
    }

    else
        return (EB_FALSE);
}

EbBool check_mv_for_tilt(uint32_t hierarchical_levels, uint32_t temporal_layer_index,
                         int32_t *x_current_mv, int32_t *y_current_mv, int32_t *x_candidate_mv,
                         int32_t *y_candidate_mv) {
    if (*x_current_mv < LOW_AMPLITUDE_TH &&
        *x_candidate_mv<
            LOW_AMPLITUDE_TH && * y_current_mv * *
            y_candidate_mv > 0 // both negative or both positives and both different than 0 i.e. same direction and non Stationary)
        && ABS(*y_current_mv) >=
               global_motion_threshold[hierarchical_levels][temporal_layer_index] // high amplitude
        && ABS(*y_candidate_mv) >=
               global_motion_threshold[hierarchical_levels][temporal_layer_index] // high amplitude
        && ABS(*y_current_mv - *y_candidate_mv) < LOW_AMPLITUDE_TH) { // close amplitude

        return (EB_TRUE);
    }

    else
        return (EB_FALSE);
}

EbBool check_mv_for_non_uniform_motion(int32_t *x_current_mv, int32_t *y_current_mv,
                                       int32_t *x_candidate_mv, int32_t *y_candidate_mv) {
    int32_t mv_threshold = 40; //LOW_AMPLITUDE_TH + 18;
    // Either the x or the y direction is greater than threshold
    if ((ABS(*x_current_mv - *x_candidate_mv) > mv_threshold) ||
        (ABS(*y_current_mv - *y_candidate_mv) > mv_threshold))
        return (EB_TRUE);
    else
        return (EB_FALSE);
}

void check_for_non_uniform_motion_vector_field(PictureParentControlSet *pcs_ptr) {
    uint32_t sb_count;
    uint32_t pic_width_in_sb =
        (pcs_ptr->enhanced_picture_ptr->width + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;
    uint32_t sb_origin_x;
    uint32_t sb_origin_y;

    int32_t  x_current_mv                   = 0;
    int32_t  y_current_mv                   = 0;
    int32_t  x_left_mv                      = 0;
    int32_t  y_left_mv                      = 0;
    int32_t  x_top_mv                       = 0;
    int32_t  y_top_mv                       = 0;
    int32_t  x_right_mv                     = 0;
    int32_t  y_right_mv                       = 0;
    int32_t  x_bottom_mv                    = 0;
    int32_t  y_bottom_mv                    = 0;
    uint32_t count_of_non_uniform_neighbors = 0;

    for (sb_count = 0; sb_count < pcs_ptr->sb_total_count; ++sb_count) {
        count_of_non_uniform_neighbors = 0;

        sb_origin_x = (sb_count % pic_width_in_sb) * BLOCK_SIZE_64;
        sb_origin_y = (sb_count / pic_width_in_sb) * BLOCK_SIZE_64;

        if (((sb_origin_x + BLOCK_SIZE_64) <= pcs_ptr->enhanced_picture_ptr->width) &&
            ((sb_origin_y + BLOCK_SIZE_64) <= pcs_ptr->enhanced_picture_ptr->height)) {
            // Current MV
            eb_get_mv(pcs_ptr, sb_count, &x_current_mv, &y_current_mv);

            // Left MV
            if (sb_origin_x == 0) {
                x_left_mv = 0;
                y_left_mv = 0;
            } else
                eb_get_mv(pcs_ptr, sb_count - 1, &x_left_mv, &y_left_mv);
            count_of_non_uniform_neighbors += check_mv_for_non_uniform_motion(
                &x_current_mv, &y_current_mv, &x_left_mv, &y_left_mv);

            // Top MV
            if (sb_origin_y == 0) {
                x_top_mv = 0;
                y_top_mv = 0;
            } else
                eb_get_mv(pcs_ptr, sb_count - pic_width_in_sb, &x_top_mv, &y_top_mv);
            count_of_non_uniform_neighbors +=
                check_mv_for_non_uniform_motion(&x_current_mv, &y_current_mv, &x_top_mv, &y_top_mv);

            // Right MV
            if ((sb_origin_x + (BLOCK_SIZE_64 << 1)) > pcs_ptr->enhanced_picture_ptr->width) {
                x_right_mv = 0;
                y_right_mv   = 0;
            } else
                eb_get_mv(pcs_ptr, sb_count + 1, &x_right_mv, &y_right_mv);
            count_of_non_uniform_neighbors += check_mv_for_non_uniform_motion(
                &x_current_mv, &y_current_mv, &x_right_mv, &y_right_mv);

            // Bottom MV
            if ((sb_origin_y + (BLOCK_SIZE_64 << 1)) > pcs_ptr->enhanced_picture_ptr->height) {
                x_bottom_mv = 0;
                y_bottom_mv = 0;
            } else
                eb_get_mv(pcs_ptr, sb_count + pic_width_in_sb, &x_bottom_mv, &y_bottom_mv);
            count_of_non_uniform_neighbors += check_mv_for_non_uniform_motion(
                &x_current_mv, &y_current_mv, &x_bottom_mv, &y_bottom_mv);
        }
    }
}

void detect_global_motion(PictureParentControlSet *pcs_ptr) {
    //initilize global motion to be OFF for all references frames.
    memset(pcs_ptr->is_global_motion, EB_FALSE, MAX_NUM_OF_REF_PIC_LIST * REF_LIST_MAX_DEPTH);

    if (pcs_ptr->gm_level <= GM_DOWN) {
        uint32_t num_of_list_to_search =
            (pcs_ptr->slice_type == P_SLICE) ? (uint32_t)REF_LIST_0 : (uint32_t)REF_LIST_1;

        for (uint32_t list_index = REF_LIST_0; list_index <= num_of_list_to_search; ++list_index) {
            uint32_t num_of_ref_pic_to_search;
            if (pcs_ptr->is_alt_ref == EB_TRUE)
                num_of_ref_pic_to_search = 1;
            else
                num_of_ref_pic_to_search = pcs_ptr->slice_type == P_SLICE
                                               ? pcs_ptr->ref_list0_count
                                               : list_index == REF_LIST_0
                                                     ? pcs_ptr->ref_list0_count
                                                     : pcs_ptr->ref_list1_count;

            // Ref Picture Loop
            for (uint32_t ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search;
                 ++ref_pic_index) {
                pcs_ptr->is_global_motion[list_index][ref_pic_index] = EB_FALSE;
                if (pcs_ptr->global_motion_estimation[list_index][ref_pic_index].wmtype >
                    TRANSLATION)
                    pcs_ptr->is_global_motion[list_index][ref_pic_index] = EB_TRUE;
            }
        }
    } else {
        uint32_t sb_count;
        uint32_t pic_width_in_sb =
            (pcs_ptr->enhanced_picture_ptr->width + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;
        uint32_t sb_origin_x;
        uint32_t sb_origin_y;

        uint32_t total_checked_sbs = 0;
        uint32_t total_pan_sbs     = 0;

        int32_t  x_current_mv   = 0;
        int32_t  y_current_mv   = 0;
        int32_t  x_left_mv      = 0;
        int32_t  y_left_mv      = 0;
        int32_t  x_top_mv       = 0;
        int32_t  y_top_mv       = 0;
        int32_t  x_right_mv     = 0;
        int32_t  y_right_mv       = 0;
        int32_t  x_bottom_mv    = 0;
        int32_t  y_bottom_mv    = 0;
        int64_t  x_tile_mv_sum  = 0;
        int64_t  y_tilt_mv_sum  = 0;
        int64_t  x_pan_mv_sum   = 0;
        int64_t  y_pan_mv_sum   = 0;
        uint32_t total_tilt_sbs = 0;

        uint32_t total_tilt_high_amp_sbs = 0;
        uint32_t total_pan_high_amp_sbs  = 0;

        for (sb_count = 0; sb_count < pcs_ptr->sb_total_count; ++sb_count) {
            sb_origin_x = (sb_count % pic_width_in_sb) * BLOCK_SIZE_64;
            sb_origin_y = (sb_count / pic_width_in_sb) * BLOCK_SIZE_64;
            if (((sb_origin_x + BLOCK_SIZE_64) <= pcs_ptr->enhanced_picture_ptr->width) &&
                ((sb_origin_y + BLOCK_SIZE_64) <= pcs_ptr->enhanced_picture_ptr->height)) {
                // Current MV
                eb_get_mv(pcs_ptr, sb_count, &x_current_mv, &y_current_mv);

                // Left MV
                if (sb_origin_x == 0) {
                    x_left_mv = 0;
                    y_left_mv = 0;
                } else
                    eb_get_mv(pcs_ptr, sb_count - 1, &x_left_mv, &y_left_mv);
                // Top MV
                if (sb_origin_y == 0) {
                    x_top_mv = 0;
                    y_top_mv = 0;
                } else
                    eb_get_mv(pcs_ptr, sb_count - pic_width_in_sb, &x_top_mv, &y_top_mv);
                // Right MV
                if ((sb_origin_x + (BLOCK_SIZE_64 << 1)) > pcs_ptr->enhanced_picture_ptr->width) {
                    x_right_mv = 0;
                    y_right_mv   = 0;
                } else
                    eb_get_mv(pcs_ptr, sb_count + 1, &x_right_mv, &y_right_mv);
                // Bottom MV
                if ((sb_origin_y + (BLOCK_SIZE_64 << 1)) > pcs_ptr->enhanced_picture_ptr->height) {
                    x_bottom_mv = 0;
                    y_bottom_mv = 0;
                } else
                    eb_get_mv(pcs_ptr, sb_count + pic_width_in_sb, &x_bottom_mv, &y_bottom_mv);
                total_checked_sbs++;

                if ((EbBool)(check_mv_for_pan(pcs_ptr->hierarchical_levels,
                                              pcs_ptr->temporal_layer_index,
                                              &x_current_mv,
                                              &y_current_mv,
                                              &x_left_mv,
                                              &y_left_mv) ||
                             check_mv_for_pan(pcs_ptr->hierarchical_levels,
                                              pcs_ptr->temporal_layer_index,
                                              &x_current_mv,
                                              &y_current_mv,
                                              &x_top_mv,
                                              &y_top_mv) ||
                             check_mv_for_pan(pcs_ptr->hierarchical_levels,
                                              pcs_ptr->temporal_layer_index,
                                              &x_current_mv,
                                              &y_current_mv,
                                              &x_right_mv,
                                              &y_right_mv) ||
                             check_mv_for_pan(pcs_ptr->hierarchical_levels,
                                              pcs_ptr->temporal_layer_index,
                                              &x_current_mv,
                                              &y_current_mv,
                                              &x_bottom_mv,
                                              &y_bottom_mv))) {
                    total_pan_sbs++;

                    x_pan_mv_sum += x_current_mv;
                    y_pan_mv_sum += y_current_mv;
                }

                if ((EbBool)(check_mv_for_tilt(pcs_ptr->hierarchical_levels,
                                               pcs_ptr->temporal_layer_index,
                                               &x_current_mv,
                                               &y_current_mv,
                                               &x_left_mv,
                                               &y_left_mv) ||
                             check_mv_for_tilt(pcs_ptr->hierarchical_levels,
                                               pcs_ptr->temporal_layer_index,
                                               &x_current_mv,
                                               &y_current_mv,
                                               &x_top_mv,
                                               &y_top_mv) ||
                             check_mv_for_tilt(pcs_ptr->hierarchical_levels,
                                               pcs_ptr->temporal_layer_index,
                                               &x_current_mv,
                                               &y_current_mv,
                                               &x_right_mv,
                                               &y_right_mv) ||
                             check_mv_for_tilt(pcs_ptr->hierarchical_levels,
                                               pcs_ptr->temporal_layer_index,
                                               &x_current_mv,
                                               &y_current_mv,
                                               &x_bottom_mv,
                                               &y_bottom_mv))) {
                    total_tilt_sbs++;

                    x_tile_mv_sum += x_current_mv;
                    y_tilt_mv_sum += y_current_mv;
                }

                if ((EbBool)(check_mv_for_pan_high_amp(pcs_ptr->hierarchical_levels,
                                                       pcs_ptr->temporal_layer_index,
                                                       &x_current_mv,
                                                       &x_left_mv) ||
                             check_mv_for_pan_high_amp(pcs_ptr->hierarchical_levels,
                                                       pcs_ptr->temporal_layer_index,
                                                       &x_current_mv,
                                                       &x_top_mv) ||
                             check_mv_for_pan_high_amp(pcs_ptr->hierarchical_levels,
                                                       pcs_ptr->temporal_layer_index,
                                                       &x_current_mv,
                                                       &x_right_mv) ||
                             check_mv_for_pan_high_amp(pcs_ptr->hierarchical_levels,
                                                       pcs_ptr->temporal_layer_index,
                                                       &x_current_mv,
                                                       &x_bottom_mv))) {
                    total_pan_high_amp_sbs++;
                }

                if ((EbBool)(check_mv_for_tilt_high_amp(pcs_ptr->hierarchical_levels,
                                                        pcs_ptr->temporal_layer_index,
                                                        &y_current_mv,
                                                        &y_left_mv) ||
                             check_mv_for_tilt_high_amp(pcs_ptr->hierarchical_levels,
                                                        pcs_ptr->temporal_layer_index,
                                                        &y_current_mv,
                                                        &y_top_mv) ||
                             check_mv_for_tilt_high_amp(pcs_ptr->hierarchical_levels,
                                                        pcs_ptr->temporal_layer_index,
                                                        &y_current_mv,
                                                        &y_right_mv) ||
                             check_mv_for_tilt_high_amp(pcs_ptr->hierarchical_levels,
                                                        pcs_ptr->temporal_layer_index,
                                                        &y_current_mv,
                                                        &y_bottom_mv))) {
                    total_tilt_high_amp_sbs++;
                }
            }
        }
        pcs_ptr->is_pan  = EB_FALSE;
        pcs_ptr->is_tilt = EB_FALSE;

        pcs_ptr->pan_mvx  = 0;
        pcs_ptr->pan_mvy  = 0;
        pcs_ptr->tilt_mvx = 0;
        pcs_ptr->tilt_mvy = 0;

        // If more than PAN_SB_PERCENTAGE % of SBs are PAN
        if ((total_checked_sbs > 0) &&
            ((total_pan_sbs * 100 / total_checked_sbs) > PAN_SB_PERCENTAGE)) {
            pcs_ptr->is_pan = EB_TRUE;
            if (total_pan_sbs > 0) {
                pcs_ptr->pan_mvx = (int16_t)(x_pan_mv_sum / total_pan_sbs);
                pcs_ptr->pan_mvy = (int16_t)(y_pan_mv_sum / total_pan_sbs);
            }
        }

        if ((total_checked_sbs > 0) &&
            ((total_tilt_sbs * 100 / total_checked_sbs) > PAN_SB_PERCENTAGE)) {
            pcs_ptr->is_tilt = EB_TRUE;
            if (total_tilt_sbs > 0) {
                pcs_ptr->tilt_mvx = (int16_t)(x_tile_mv_sum / total_tilt_sbs);
                pcs_ptr->tilt_mvy = (int16_t)(y_tilt_mv_sum / total_tilt_sbs);
            }
        }
    }
}

static void initial_rate_control_context_dctor(EbPtr p) {
    EbThreadContext *          thread_context_ptr = (EbThreadContext *)p;
    InitialRateControlContext *obj = (InitialRateControlContext *)thread_context_ptr->priv;
    EB_FREE_ARRAY(obj);
}

/************************************************
* Initial Rate Control Context Constructor
************************************************/
EbErrorType initial_rate_control_context_ctor(EbThreadContext *  thread_context_ptr,
                                              const EbEncHandle *enc_handle_ptr) {
    InitialRateControlContext *context_ptr;
    EB_CALLOC_ARRAY(context_ptr, 1);
    thread_context_ptr->priv  = context_ptr;
    thread_context_ptr->dctor = initial_rate_control_context_dctor;

    context_ptr->motion_estimation_results_input_fifo_ptr = eb_system_resource_get_consumer_fifo(
        enc_handle_ptr->motion_estimation_results_resource_ptr, 0);
    context_ptr->initialrate_control_results_output_fifo_ptr = eb_system_resource_get_producer_fifo(
        enc_handle_ptr->initial_rate_control_results_resource_ptr, 0);

    return EB_ErrorNone;
}

/************************************************
* Release Pa Reference Objects
** Check if reference pictures are needed
** release them when appropriate
************************************************/
void release_pa_reference_objects(SequenceControlSet *scs_ptr, PictureParentControlSet *pcs_ptr) {
    // PA Reference Pictures
    uint32_t num_of_list_to_search;
    uint32_t list_index;
    uint32_t ref_pic_index;
    if (pcs_ptr->slice_type != I_SLICE) {
        num_of_list_to_search = (pcs_ptr->slice_type == P_SLICE) ? REF_LIST_0 : REF_LIST_1;

        // List Loop
        for (list_index = REF_LIST_0; list_index <= num_of_list_to_search; ++list_index) {
            // Release PA Reference Pictures
            uint8_t num_of_ref_pic_to_search =
                (pcs_ptr->slice_type == P_SLICE)
                    ? MIN(pcs_ptr->ref_list0_count, scs_ptr->reference_count)
                    : (list_index == REF_LIST_0)
                          ? MIN(pcs_ptr->ref_list0_count, scs_ptr->reference_count)
                          : MIN(pcs_ptr->ref_list1_count, scs_ptr->reference_count);

            for (ref_pic_index = 0; ref_pic_index < num_of_ref_pic_to_search; ++ref_pic_index) {
                if (pcs_ptr->ref_pa_pic_ptr_array[list_index][ref_pic_index] != EB_NULL) {
                    eb_release_object(pcs_ptr->ref_pa_pic_ptr_array[list_index][ref_pic_index]);
                }
            }
        }
    }

    if (pcs_ptr->pa_reference_picture_wrapper_ptr != EB_NULL) {
        eb_release_object(pcs_ptr->pa_reference_picture_wrapper_ptr);
    }

    return;
}

/************************************************
* Global Motion Detection Based on ME information
** Mark pictures for pan
** Mark pictures for tilt
** No lookahead information used in this function
************************************************/
void me_based_global_motion_detection(PictureParentControlSet *pcs_ptr) {
    // PAN Generation
    pcs_ptr->is_pan  = EB_FALSE;
    pcs_ptr->is_tilt = EB_FALSE;

    if (pcs_ptr->slice_type != I_SLICE) detect_global_motion(pcs_ptr);
    // Check if the motion vector field for temporal layer 0 pictures
    if (pcs_ptr->slice_type != I_SLICE && pcs_ptr->temporal_layer_index == 0)
        check_for_non_uniform_motion_vector_field(pcs_ptr);
    return;
}
/************************************************
* Global Motion Detection Based on Lookahead
** Mark pictures for pan
** Mark pictures for tilt
** LAD Window: min (8 or sliding window size)
************************************************/
void update_global_motion_detection_over_time(EncodeContext *          encode_context_ptr,
                                              SequenceControlSet *     scs_ptr,
                                              PictureParentControlSet *pcs_ptr) {
    InitialRateControlReorderEntry *temp_queue_entry_ptr;
    PictureParentControlSet *       temp_pcs_ptr;

    uint32_t total_pan_pictures     = 0;
    uint32_t total_checked_pictures = 0;
    uint32_t total_tilt_pictures    = 0;
    uint32_t update_is_pan_frames_to_check;
    uint32_t input_queue_index;
    uint32_t frames_to_check_index;

    (void)scs_ptr;

    // Determine number of frames to check (8 frames)
    update_is_pan_frames_to_check = MIN(8, pcs_ptr->frames_in_sw);

    // Walk the first N entries in the sliding window
    input_queue_index = encode_context_ptr->initial_rate_control_reorder_queue_head_index;
    uint32_t update_frames_to_check = update_is_pan_frames_to_check;
    for (frames_to_check_index = 0; frames_to_check_index < update_frames_to_check;
         frames_to_check_index++) {
        temp_queue_entry_ptr =
            encode_context_ptr->initial_rate_control_reorder_queue[input_queue_index];
        temp_pcs_ptr =
            ((PictureParentControlSet *)(temp_queue_entry_ptr->parent_pcs_wrapper_ptr)->object_ptr);

        if (temp_pcs_ptr->slice_type != I_SLICE) {
            total_pan_pictures += (temp_pcs_ptr->is_pan == EB_TRUE);

            total_tilt_pictures += (temp_pcs_ptr->is_tilt == EB_TRUE);

            // Keep track of checked pictures
            total_checked_pictures++;
        }

        // Increment the input_queue_index Iterator
        input_queue_index = (input_queue_index == INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                                ? 0
                                : input_queue_index + 1;
    }

    pcs_ptr->is_pan  = EB_FALSE;
    pcs_ptr->is_tilt = EB_FALSE;

    if (total_checked_pictures) {
        if (pcs_ptr->slice_type != I_SLICE) {
            if ((total_pan_pictures * 100 / total_checked_pictures) > 75) pcs_ptr->is_pan = EB_TRUE;
        }
    }
    return;
}

/************************************************
* Update BEA Information Based on Lookahead
** Average zzCost of Collocated SB throughout lookahead frames
** Set isMostOfPictureNonMoving based on number of non moving SBs
** LAD Window: min (2xmgpos+1 or sliding window size)
************************************************/

void update_bea_info_over_time(EncodeContext *          encode_context_ptr,
                               PictureParentControlSet *pcs_ptr) {
    InitialRateControlReorderEntry *temp_queue_entry_ptr;
    PictureParentControlSet *       temp_pcs_ptr;
    uint32_t                        update_non_moving_index_array_frames_to_check;
    uint16_t                        sb_idx;
    uint16_t                        frames_to_check_index;
    uint64_t                        non_moving_index_sum = 0;
    uint32_t                        input_queue_index;

    SequenceControlSet *scs_ptr = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
    // Update motionIndexArray of the current picture by averaging the motionIndexArray of the N future pictures
    // Determine number of frames to check N
    update_non_moving_index_array_frames_to_check =
        MIN(MIN(((pcs_ptr->pred_struct_ptr->pred_struct_period << 1) + 1), pcs_ptr->frames_in_sw),
            scs_ptr->static_config.look_ahead_distance);
    uint64_t me_dist           = 0;
    uint8_t  me_dist_pic_count = 0;
    // SB Loop
    for (sb_idx = 0; sb_idx < pcs_ptr->sb_total_count; ++sb_idx) {
        uint16_t non_moving_index_over_sliding_window = pcs_ptr->non_moving_index_array[sb_idx];

        // Walk the first N entries in the sliding window starting picture + 1
        input_queue_index =
            (encode_context_ptr->initial_rate_control_reorder_queue_head_index ==
             INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                ? 0
                : encode_context_ptr->initial_rate_control_reorder_queue_head_index + 1;
        for (frames_to_check_index = 0;
             frames_to_check_index < update_non_moving_index_array_frames_to_check - 1;
             frames_to_check_index++) {
            temp_queue_entry_ptr =
                encode_context_ptr->initial_rate_control_reorder_queue[input_queue_index];
            temp_pcs_ptr =
                ((PictureParentControlSet *)(temp_queue_entry_ptr->parent_pcs_wrapper_ptr)
                     ->object_ptr);

            if (temp_pcs_ptr->slice_type == I_SLICE || temp_pcs_ptr->end_of_sequence_flag) break;
            // Limit the distortion to lower layers 0, 1 and 2 only. Higher layers have close temporal distance and lower distortion that might contaminate the data
            if (temp_pcs_ptr->temporal_layer_index <
                MAX((int8_t)pcs_ptr->hierarchical_levels - 1, 2)) {
                if (sb_idx == 0) me_dist_pic_count++;
                me_dist += (temp_pcs_ptr->slice_type == I_SLICE)
                               ? 0
                               : (uint64_t)temp_pcs_ptr->rc_me_distortion[sb_idx];
            }
            // Store the filtered_sse of next ALT_REF picture in the I slice to be used in QP Scaling
            if (pcs_ptr->slice_type == I_SLICE && pcs_ptr->filtered_sse == 0 && sb_idx == 0 &&
                temp_pcs_ptr->temporal_layer_index == 0) {
                pcs_ptr->filtered_sse    = temp_pcs_ptr->filtered_sse;
                pcs_ptr->filtered_sse_uv = temp_pcs_ptr->filtered_sse_uv;
            }
            non_moving_index_over_sliding_window += temp_pcs_ptr->non_moving_index_array[sb_idx];

            // Increment the input_queue_index Iterator
            input_queue_index =
                (input_queue_index == INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                    ? 0
                    : input_queue_index + 1;
        }
        pcs_ptr->non_moving_index_array[sb_idx] =
            (uint8_t)(non_moving_index_over_sliding_window / (frames_to_check_index + 1));

        non_moving_index_sum += pcs_ptr->non_moving_index_array[sb_idx];
    }

    pcs_ptr->non_moving_index_average = (uint16_t)non_moving_index_sum / pcs_ptr->sb_total_count;
    me_dist_pic_count                 = MAX(me_dist_pic_count, 1);
    pcs_ptr->qp_scaling_average_complexity =
        (uint16_t)((uint64_t)me_dist / pcs_ptr->sb_total_count / 256 / me_dist_pic_count);
    return;
}

/****************************************
* Init ZZ Cost array to default values
** Used when no Lookahead is available
****************************************/
void init_zz_cost_info(PictureParentControlSet *pcs_ptr) {
    uint16_t sb_idx;
    pcs_ptr->non_moving_index_average = INVALID_ZZ_COST;

    // SB Loop
    for (sb_idx = 0; sb_idx < pcs_ptr->sb_total_count; ++sb_idx)
        pcs_ptr->non_moving_index_array[sb_idx] = INVALID_ZZ_COST;
    return;
}

/************************************************
* Update uniform motion field
** Update Uniformly moving SBs using
** collocated SBs infor in lookahead pictures
** LAD Window: min (2xmgpos+1 or sliding window size)
************************************************/
void update_motion_field_uniformity_over_time(EncodeContext *          encode_context_ptr,
                                              SequenceControlSet *     scs_ptr,
                                              PictureParentControlSet *pcs_ptr) {
    InitialRateControlReorderEntry *temp_queue_entry_ptr;
    PictureParentControlSet *       temp_pcs_ptr;
    uint32_t                        input_queue_index;
    uint32_t                        no_frames_to_check;
    uint32_t                        frames_to_check_index;
    //SVT_LOG("To update POC %d\tframesInSw = %d\n", pcs_ptr->picture_number, pcs_ptr->frames_in_sw);
    // Determine number of frames to check N
    no_frames_to_check =
        MIN(MIN(((pcs_ptr->pred_struct_ptr->pred_struct_period << 1) + 1), pcs_ptr->frames_in_sw),
            scs_ptr->static_config.look_ahead_distance);

    // Walk the first N entries in the sliding window starting picture + 1
    input_queue_index = (encode_context_ptr->initial_rate_control_reorder_queue_head_index ==
                         INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                            ? 0
                            : encode_context_ptr->initial_rate_control_reorder_queue_head_index;
    for (frames_to_check_index = 0; frames_to_check_index < no_frames_to_check - 1;
         frames_to_check_index++) {
        temp_queue_entry_ptr =
            encode_context_ptr->initial_rate_control_reorder_queue[input_queue_index];
        temp_pcs_ptr =
            ((PictureParentControlSet *)(temp_queue_entry_ptr->parent_pcs_wrapper_ptr)->object_ptr);

        if (temp_pcs_ptr->end_of_sequence_flag) break;
        // Increment the input_queue_index Iterator
        input_queue_index = (input_queue_index == INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                                ? 0
                                : input_queue_index + 1;
    }
    return;
}
InitialRateControlReorderEntry *determine_picture_offset_in_queue(
    EncodeContext *encode_context_ptr, PictureParentControlSet *pcs_ptr,
    MotionEstimationResults *in_results_ptr) {
    InitialRateControlReorderEntry *queue_entry_ptr;
    int32_t                         queue_entry_index;

    queue_entry_index =
        (int32_t)(pcs_ptr->picture_number -
                  encode_context_ptr
                      ->initial_rate_control_reorder_queue
                          [encode_context_ptr->initial_rate_control_reorder_queue_head_index]
                      ->picture_number);
    queue_entry_index += encode_context_ptr->initial_rate_control_reorder_queue_head_index;
    queue_entry_index = (queue_entry_index > INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                            ? queue_entry_index - INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH
                            : queue_entry_index;
    queue_entry_ptr = encode_context_ptr->initial_rate_control_reorder_queue[queue_entry_index];
    queue_entry_ptr->parent_pcs_wrapper_ptr = in_results_ptr->pcs_wrapper_ptr;
    queue_entry_ptr->picture_number         = pcs_ptr->picture_number;

    return queue_entry_ptr;
}

void get_histogram_queue_data(SequenceControlSet *scs_ptr, EncodeContext *encode_context_ptr,
                              PictureParentControlSet *pcs_ptr) {
    HlRateControlHistogramEntry *histogram_queue_entry_ptr;
    int32_t                      histogram_queue_entry_index;

    // Determine offset from the Head Ptr for HLRC histogram queue
    eb_block_on_mutex(scs_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
    histogram_queue_entry_index = (int32_t)(
        pcs_ptr->picture_number -
        encode_context_ptr
            ->hl_rate_control_historgram_queue[encode_context_ptr
                                                   ->hl_rate_control_historgram_queue_head_index]
            ->picture_number);
    histogram_queue_entry_index += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
    histogram_queue_entry_index =
        (histogram_queue_entry_index > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1)
            ? histogram_queue_entry_index - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH
            : (histogram_queue_entry_index < 0)
                  ? histogram_queue_entry_index + HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH
                  : histogram_queue_entry_index;
    histogram_queue_entry_ptr =
        encode_context_ptr->hl_rate_control_historgram_queue[histogram_queue_entry_index];

    //histogram_queue_entry_ptr->parent_pcs_wrapper_ptr  = in_results_ptr->pcs_wrapper_ptr;
    histogram_queue_entry_ptr->picture_number       = pcs_ptr->picture_number;
    histogram_queue_entry_ptr->end_of_sequence_flag = pcs_ptr->end_of_sequence_flag;
    histogram_queue_entry_ptr->slice_type           = pcs_ptr->slice_type;
    histogram_queue_entry_ptr->temporal_layer_index = pcs_ptr->temporal_layer_index;
    histogram_queue_entry_ptr->full_sb_count        = pcs_ptr->full_sb_count;
    histogram_queue_entry_ptr->life_count           = 0;
    histogram_queue_entry_ptr->passed_to_hlrc       = EB_FALSE;
    histogram_queue_entry_ptr->is_coded             = EB_FALSE;
    histogram_queue_entry_ptr->total_num_bits_coded = 0;
    histogram_queue_entry_ptr->frames_in_sw         = 0;
    EB_MEMCPY(histogram_queue_entry_ptr->me_distortion_histogram,
              pcs_ptr->me_distortion_histogram,
              sizeof(uint16_t) * NUMBER_OF_SAD_INTERVALS);

    EB_MEMCPY(histogram_queue_entry_ptr->ois_distortion_histogram,
              pcs_ptr->ois_distortion_histogram,
              sizeof(uint16_t) * NUMBER_OF_INTRA_SAD_INTERVALS);

    eb_release_mutex(scs_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);
    //SVT_LOG("Test1 POC: %d\t POC: %d\t LifeCount: %d\n", histogram_queue_entry_ptr->picture_number, pcs_ptr->picture_number,  histogram_queue_entry_ptr->life_count);

    return;
}

void update_histogram_queue_entry(SequenceControlSet *scs_ptr, EncodeContext *encode_context_ptr,
                                  PictureParentControlSet *pcs_ptr, uint32_t frames_in_sw) {
    HlRateControlHistogramEntry *histogram_queue_entry_ptr;
    int32_t                      histogram_queue_entry_index;

    eb_block_on_mutex(scs_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);

    histogram_queue_entry_index = (int32_t)(
        pcs_ptr->picture_number -
        encode_context_ptr
            ->hl_rate_control_historgram_queue[encode_context_ptr
                                                   ->hl_rate_control_historgram_queue_head_index]
            ->picture_number);
    histogram_queue_entry_index += encode_context_ptr->hl_rate_control_historgram_queue_head_index;
    histogram_queue_entry_index =
        (histogram_queue_entry_index > HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH - 1)
            ? histogram_queue_entry_index - HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH
            : (histogram_queue_entry_index < 0)
                  ? histogram_queue_entry_index + HIGH_LEVEL_RATE_CONTROL_HISTOGRAM_QUEUE_MAX_DEPTH
                  : histogram_queue_entry_index;
    histogram_queue_entry_ptr =
        encode_context_ptr->hl_rate_control_historgram_queue[histogram_queue_entry_index];
    histogram_queue_entry_ptr->passed_to_hlrc = EB_TRUE;

    if (scs_ptr->static_config.rate_control_mode == 2)
        histogram_queue_entry_ptr->life_count +=
            (int16_t)(scs_ptr->static_config.intra_period_length + 1) -
            3; // FramelevelRC does not decrease the life count for first picture in each temporal layer
    else
        histogram_queue_entry_ptr->life_count += pcs_ptr->historgram_life_count;

    histogram_queue_entry_ptr->frames_in_sw = frames_in_sw;
    eb_release_mutex(scs_ptr->encode_context_ptr->hl_rate_control_historgram_queue_mutex);

    return;
}

// Derives blockinessPresentFlag
void DeriveBlockinessPresentFlag(
    SequenceControlSet        *sequence_control_set_ptr,
    PictureParentControlSet   *picture_control_set_ptr)
{
    uint32_t                      sb_index;

    for (sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
        SbParams         *lcuParamPtr = &sequence_control_set_ptr->sb_params_array[sb_index];
        picture_control_set_ptr->complex_sb_array[sb_index] = SB_COMPLEXITY_STATUS_INVALID;

        // Spatially complex SB within a spatially complex area
        if (IsSpatiallyComplexArea(
            picture_control_set_ptr,
            (picture_control_set_ptr->enhanced_picture_ptr->width + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64,
            sb_index,
            lcuParamPtr->origin_x,
            lcuParamPtr->origin_y)) {
            // Active SB within an active scene (added a check on 4K & non-BASE to restrict the action - could be generated for all resolutions/layers)
            if (picture_control_set_ptr->non_moving_index_array[sb_index] == SB_COMPLEXITY_NON_MOVING_INDEX_TH_0 && picture_control_set_ptr->non_moving_index_average >= SB_COMPLEXITY_NON_MOVING_INDEX_TH_1 && picture_control_set_ptr->temporal_layer_index > 0 && sequence_control_set_ptr->input_resolution == INPUT_SIZE_4K_RANGE)
                picture_control_set_ptr->complex_sb_array[sb_index] = SB_COMPLEXITY_STATUS_2;
            // Active SB within a scene with a moderate acitivity (eg. active foregroud & static background)
            else if (picture_control_set_ptr->non_moving_index_array[sb_index] == SB_COMPLEXITY_NON_MOVING_INDEX_TH_0 && picture_control_set_ptr->non_moving_index_average >= SB_COMPLEXITY_NON_MOVING_INDEX_TH_2 && picture_control_set_ptr->non_moving_index_average < SB_COMPLEXITY_NON_MOVING_INDEX_TH_1)
                picture_control_set_ptr->complex_sb_array[sb_index] = SB_COMPLEXITY_STATUS_1;
            else
                picture_control_set_ptr->complex_sb_array[sb_index] = SB_COMPLEXITY_STATUS_0;
        }
        else
            picture_control_set_ptr->complex_sb_array[sb_index] = SB_COMPLEXITY_STATUS_0;
    }
}

#if CUTREE_LA
static AOM_INLINE void get_quantize_error(MacroblockPlane *p, int plane,
                                          tran_low_t *coeff, tran_low_t *qcoeff,
                                          tran_low_t *dqcoeff, TxSize tx_size,
                                          uint16_t *eob, int64_t *recon_error,
                                          int64_t *sse) {
  const ScanOrder *const scan_order = &av1_scan_orders[tx_size][DCT_1D]; //&av1_default_scan_orders[tx_size]
  int pix_num = 1 << num_pels_log2_lookup[txsize_to_bsize[tx_size]];
  const int shift = tx_size == TX_32X32 ? 0 : 2;

  eb_av1_quantize_fp(coeff, pix_num, p->zbin_QTX, p->round_fp_QTX, p->quant_fp_QTX,
                  p->quant_shift_QTX, qcoeff, dqcoeff, p->dequant_QTX, eob,
                  scan_order->scan, scan_order->iscan);

  *recon_error = av1_block_error(coeff, dqcoeff, pix_num, sse) >> shift;
  *recon_error = AOMMAX(*recon_error, 1);

  *sse = (*sse) >> shift;
  *sse = AOMMAX(*sse, 1);
}

static int rate_estimator(tran_low_t *qcoeff, int eob, TxSize tx_size) {
  const ScanOrder *const scan_order = &av1_scan_orders[tx_size][DCT_1D]; //&av1_default_scan_orders[tx_size]

  assert((1 << num_pels_log2_lookup[txsize_to_bsize[tx_size]]) >= eob);

  int rate_cost = 1;

  for (int idx = 0; idx < eob; ++idx) {
    int abs_level = abs(qcoeff[scan_order->scan[idx]]);
    rate_cost += (int)(log(abs_level + 1.0) / log(2.0)) + 1;
  }

  return (rate_cost << AV1_PROB_COST_SHIFT);
}

static void result_model_store(PictureParentControlSet *picture_control_set_ptr, OisMbResults *ois_mb_results_ptr,
        uint32_t mb_origin_x, uint32_t mb_origin_y, uint32_t picture_width_in_mb) {
  const int mi_height = mi_size_high[TX_16X16];
  const int mi_width = mi_size_wide[TX_16X16];
  const int step = 1 << 2;//cpi->tpl_stats_block_mis_log2; //is_720p_or_larger ? 2 : 1;

  int64_t intra_cost = ois_mb_results_ptr->intra_cost / (mi_height * mi_width);
  int64_t inter_cost = ois_mb_results_ptr->inter_cost / (mi_height * mi_width);
  int64_t srcrf_dist = ois_mb_results_ptr->srcrf_dist / (mi_height * mi_width);
  int64_t recrf_dist = ois_mb_results_ptr->recrf_dist / (mi_height * mi_width);
  int64_t srcrf_rate = ois_mb_results_ptr->srcrf_rate / (mi_height * mi_width);
  int64_t recrf_rate = ois_mb_results_ptr->recrf_rate / (mi_height * mi_width);

  intra_cost = AOMMAX(1, intra_cost);
  inter_cost = AOMMAX(1, inter_cost);
  srcrf_dist = AOMMAX(1, srcrf_dist);
  recrf_dist = AOMMAX(1, recrf_dist);
  srcrf_rate = AOMMAX(1, srcrf_rate);
  recrf_rate = AOMMAX(1, recrf_rate);

  for (int idy = 0; idy < mi_height; idy += step) {
    //TplDepStats *tpl_ptr =
    //    &tpl_stats_ptr[av1_tpl_ptr_pos(cpi, mi_row + idy, mi_col, stride)];
    OisMbResults *dst_ptr = picture_control_set_ptr->ois_mb_results[(mb_origin_y >> 4) * picture_width_in_mb + (mb_origin_x >> 4)];
    for (int idx = 0; idx < mi_width; idx += step) {
      dst_ptr->intra_cost = intra_cost;
      dst_ptr->inter_cost = inter_cost;
      dst_ptr->srcrf_dist = srcrf_dist;
      dst_ptr->recrf_dist = recrf_dist;
      dst_ptr->srcrf_rate = srcrf_rate;
      dst_ptr->recrf_rate = recrf_rate;
      dst_ptr->mv = ois_mb_results_ptr->mv;
      dst_ptr->ref_frame_poc = ois_mb_results_ptr->ref_frame_poc;
      ++dst_ptr;
    }
  }
}

#define TPL_DEP_COST_SCALE_LOG2 4
/************************************************
* Genrate CUTree MC Flow Dispenser  Based on Lookahead
** LAD Window: sliding window size
************************************************/
void cutree_mc_flow_dispenser(
    EncodeContext                   *encode_context_ptr,
    SequenceControlSet              *sequence_control_set_ptr,
    PictureParentControlSet         *picture_control_set_ptr)
{
    InitialRateControlReorderEntry   *temporaryQueueEntryPtr;
    PictureParentControlSet          *temporaryPictureControlSetPtr;

    uint32_t    inputQueueIndex;
    uint32_t    frame_to_check_index;
    uint32_t    picture_width_in_sb = (picture_control_set_ptr->enhanced_picture_ptr->width + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;
    uint32_t    picture_width_in_mb = (picture_control_set_ptr->enhanced_picture_ptr->width + 16 - 1) / 16;
    uint32_t    picture_height_in_sb = (picture_control_set_ptr->enhanced_picture_ptr->height + BLOCK_SIZE_64 - 1) / BLOCK_SIZE_64;
    uint32_t    sb_origin_x;
    uint32_t    sb_origin_y;
    int32_t     x_curr_mv = 0;
    int32_t     y_curr_mv = 0;
    uint32_t    me_mb_offset = 0;
    TxSize      tx_size = TX_16X16;
    MeCandidate *me_block_result;
    EbPictureBufferDesc  *ref_pic_ptr;
    EbReferenceObject    *referenceObject;
    struct      ScaleFactors sf;
    BlockGeom   blk_geom;
    Av1Common *cm = picture_control_set_ptr->av1_cm;
    uint32_t    kernel = (EIGHTTAP_REGULAR << 16) | EIGHTTAP_REGULAR;
    EbPictureBufferDesc *input_picture_ptr = picture_control_set_ptr->enhanced_picture_ptr;
    int64_t recon_error = 1, sse = 1;

    (void)sequence_control_set_ptr;

    DECLARE_ALIGNED(32, uint8_t, predictor8[256 * 2]);
    DECLARE_ALIGNED(32, int16_t, src_diff[256]);
    DECLARE_ALIGNED(32, tran_low_t, coeff[256]);
    DECLARE_ALIGNED(32, tran_low_t, qcoeff[256]);
    DECLARE_ALIGNED(32, tran_low_t, dqcoeff[256]);
    DECLARE_ALIGNED(32, tran_low_t, best_coeff[256]);
    uint8_t *predictor = predictor8;

    blk_geom.bwidth  = 16;
    blk_geom.bheight = 16;
    blk_geom.origin_x = 0;
    blk_geom.origin_y = 0;


    av1_setup_scale_factors_for_frame(
                &sf, picture_width_in_sb * BLOCK_SIZE_64,
                picture_height_in_sb * BLOCK_SIZE_64,
                picture_width_in_sb * BLOCK_SIZE_64,
                picture_height_in_sb * BLOCK_SIZE_64);

    MacroblockPlane mb_plane;
    int32_t qIndex = quantizer_to_qindex[(uint8_t)sequence_control_set_ptr->qp];
    Quants *const quantsMd = &picture_control_set_ptr->quantsMd;
    Dequants *const dequantsMd = &picture_control_set_ptr->deqMd;
    eb_av1_set_quantizer(
        picture_control_set_ptr,
        picture_control_set_ptr->frm_hdr.quantization_params.base_q_idx);
    eb_av1_build_quantizer(
        /*picture_control_set_ptr->hbd_mode_decision ? AOM_BITS_10 :*/ AOM_BITS_8,
        picture_control_set_ptr->frm_hdr.quantization_params.delta_q_dc[AOM_PLANE_Y],
        picture_control_set_ptr->frm_hdr.quantization_params.delta_q_dc[AOM_PLANE_U],
        picture_control_set_ptr->frm_hdr.quantization_params.delta_q_ac[AOM_PLANE_U],
        picture_control_set_ptr->frm_hdr.quantization_params.delta_q_dc[AOM_PLANE_V],
        picture_control_set_ptr->frm_hdr.quantization_params.delta_q_ac[AOM_PLANE_V],
        quantsMd,
        dequantsMd);
    mb_plane.quant_QTX       = picture_control_set_ptr->quantsMd.y_quant[qIndex];
    mb_plane.quant_fp_QTX    = picture_control_set_ptr->quantsMd.y_quant_fp[qIndex];
    mb_plane.round_fp_QTX    = picture_control_set_ptr->quantsMd.y_round_fp[qIndex];
    mb_plane.quant_shift_QTX = picture_control_set_ptr->quantsMd.y_quant_shift[qIndex];
    mb_plane.zbin_QTX        = picture_control_set_ptr->quantsMd.y_zbin[qIndex];
    mb_plane.round_QTX       = picture_control_set_ptr->quantsMd.y_round[qIndex];
    mb_plane.dequant_QTX     = picture_control_set_ptr->deqMd.y_dequant_QTX[qIndex];

    // Walk the first N entries in the sliding window
    inputQueueIndex = encode_context_ptr->initial_rate_control_reorder_queue_head_index;
    //inputQueueIndex = (inputQueueIndex == INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1) ? 0 : inputQueueIndex + 1;
    for (uint32_t sb_index = 0; sb_index < picture_control_set_ptr->sb_total_count; ++sb_index) {
        sb_origin_x = (sb_index % picture_width_in_sb) * BLOCK_SIZE_64;
        sb_origin_y = (sb_index / picture_width_in_sb) * BLOCK_SIZE_64;
        if (((sb_origin_x + BLOCK_SIZE_64) <= input_picture_ptr->width) &&
            ((sb_origin_y + BLOCK_SIZE_64) <= input_picture_ptr->height)) {
            SbParams *sb_params = &sequence_control_set_ptr->sb_params_array[sb_index];
            uint32_t pa_blk_index = 0;
            while (pa_blk_index < CU_MAX_COUNT) {
                const CodedUnitStats *blk_stats_ptr;
                blk_stats_ptr = get_coded_unit_stats(pa_blk_index);
                uint8_t bsize = blk_stats_ptr->size;
                if(bsize != 16) {
                    pa_blk_index++;
                    continue;
                }
                if (sb_params->raster_scan_cu_validity[md_scan_to_raster_scan[pa_blk_index]]) {
                    uint32_t mb_origin_x = sb_params->origin_x + blk_stats_ptr->origin_x;
                    uint32_t mb_origin_y = sb_params->origin_y + blk_stats_ptr->origin_y;
//printf("kelvin begin sb_index=%d, mb_origin_x=%d, mb_origin_y=%d, pa_blk_index=%d\n", sb_index, mb_origin_x, mb_origin_y, pa_blk_index);
                    int64_t inter_cost;
                    int32_t best_rf_idx = -1;
                    int64_t best_inter_cost = INT64_MAX;
                    MV final_best_mv = {0, 0};
                    uint32_t max_inter_ref = ((sequence_control_set_ptr->mrp_mode == 0) ? ME_MV_MRP_MODE_0 : ME_MV_MRP_MODE_1);
                    OisMbResults *ois_mb_results_ptr = picture_control_set_ptr->ois_mb_results[(mb_origin_y >> 4) * picture_width_in_mb + (mb_origin_x >> 4)];
                    int64_t best_intra_cost = ois_mb_results_ptr->intra_cost;
                    uint8_t best_mode = DC_PRED;
                    //uint8_t *src_mb = input_picture_ptr->buffer_y + mb_origin_y * input_picture_ptr->stride_y + mb_origin_x;
                    uint8_t *src_mb = input_picture_ptr->buffer_y + input_picture_ptr->origin_x + mb_origin_x +
                                     (input_picture_ptr->origin_y + mb_origin_y) * input_picture_ptr->stride_y;
                    for(uint32_t rf_idx = 0; rf_idx < max_inter_ref; rf_idx++) {
                        me_mb_offset = get_me_info_index(picture_control_set_ptr->max_number_of_pus_per_sb, &blk_geom, 0, 0);
                        me_block_result = GetMbMv(picture_control_set_ptr, rf_idx, sb_index, me_mb_offset, &x_curr_mv, &y_curr_mv);

                        uint32_t list_index = (sequence_control_set_ptr->mrp_mode == 0) ? (rf_idx < 4 ? 0 : 1)
                                                                                        : (rf_idx < 2 ? 0 : 1);
                        uint32_t ref_pic_index = (sequence_control_set_ptr->mrp_mode == 0) ? (rf_idx >= 4 ? (rf_idx - 4) : rf_idx)
                                                                                           : (rf_idx >= 2 ? (rf_idx - 2) : rf_idx);
                        if(!picture_control_set_ptr->ref_pa_pic_ptr_array[list_index][ref_pic_index])
                            continue;
                        referenceObject = (EbReferenceObject*)picture_control_set_ptr->ref_pa_pic_ptr_array[list_index][ref_pic_index]->object_ptr;
                        ref_pic_ptr = /*is16bit ? (EbPictureBufferDesc*)referenceObject->reference_picture16bit : */(EbPictureBufferDesc*)referenceObject->reference_picture;
                        const int ref_basic_offset = ref_pic_ptr->origin_y * ref_pic_ptr->stride_y + ref_pic_ptr->origin_x;
                        const int ref_mb_offset = mb_origin_y * ref_pic_ptr->stride_y + mb_origin_x;
                        uint8_t *ref_mb = ref_pic_ptr->buffer_y + ref_basic_offset + ref_mb_offset;

                        struct buf_2d ref_buf = { NULL, ref_pic_ptr->buffer_y + ref_basic_offset,
                                                  ref_pic_ptr->width, ref_pic_ptr->height,
                                                  ref_pic_ptr->stride_y };
                        InterPredParams inter_pred_params;
                        av1_init_inter_params(&inter_pred_params, 16, 16, mb_origin_y,
                                mb_origin_x, 0, 0, 8/*xd->bd*/, 0/*is_cur_buf_hbd(xd)*/, 0,
                                &sf, &ref_buf, kernel);

                        inter_pred_params.conv_params = get_conv_params(0, 0, 0, 8/*xd->bd*/);

                        MV best_mv = {x_curr_mv, y_curr_mv};
                        av1_build_inter_predictor(ref_mb, input_picture_ptr->stride_y, predictor, 16,
                                &best_mv, mb_origin_x, mb_origin_y, &inter_pred_params);
                        aom_subtract_block(16, 16, src_diff, 16, src_mb, input_picture_ptr->stride_y, predictor, 16);

                        wht_fwd_txfm(src_diff, 16, coeff, tx_size, 8/*xd->bd*/, 0/*is_cur_buf_hbd(xd)*/);

                        inter_cost = aom_satd(coeff, 256);
                        if (inter_cost < best_inter_cost) {
                            memcpy(best_coeff, coeff, sizeof(best_coeff));
                            best_rf_idx = rf_idx;
                            best_inter_cost = inter_cost;
                            final_best_mv = best_mv;

                            if (best_inter_cost < best_intra_cost) best_mode = NEWMV;
                        }
                    } // rf_idx
                    if(best_inter_cost < INT64_MAX) {
                        uint16_t eob;
                        get_quantize_error(&mb_plane, 0, best_coeff, qcoeff, dqcoeff, tx_size, &eob, &recon_error, &sse);
                        int rate_cost = rate_estimator(qcoeff, eob, tx_size);
                        ois_mb_results_ptr->srcrf_rate = rate_cost << TPL_DEP_COST_SCALE_LOG2;
                    }
                    best_intra_cost = AOMMAX(best_intra_cost, 1);
//printf("kelvincost1 poc%d sb_index=%d, mb_origin_xy=%d %d, best_mode=%d, best_intra_cost=%d, offset=%d\n", picture_control_set_ptr->picture_number, sb_index, mb_origin_x, mb_origin_y, ois_mb_results_ptr->intra_mode, best_intra_cost, (mb_origin_y >> 4) * picture_width_in_mb + (mb_origin_x >> 4));
                    if (0)//(frame_idx == 0)
                        best_inter_cost = 0;
                    else
                        best_inter_cost = AOMMIN(best_intra_cost, best_inter_cost);
                    ois_mb_results_ptr->inter_cost = best_inter_cost << TPL_DEP_COST_SCALE_LOG2;
                    ois_mb_results_ptr->intra_cost = best_intra_cost << TPL_DEP_COST_SCALE_LOG2;

                    ois_mb_results_ptr->srcrf_dist = recon_error << (TPL_DEP_COST_SCALE_LOG2);

                    if (best_mode == NEWMV) {
                        // inter recon
                        struct buf_2d ref_buf = { NULL, input_picture_ptr->buffer_y,
                            input_picture_ptr->width, input_picture_ptr->width,
                            input_picture_ptr->stride_y };
                        InterPredParams inter_pred_params;
                        av1_init_inter_params(&inter_pred_params, 16, 16, mb_origin_y,
                            mb_origin_x, 0, 0, 8/*xd->bd*/, 0/*is_cur_buf_hbd(xd)*/, 0,
                            &sf, &ref_buf, kernel);
                        inter_pred_params.conv_params = get_conv_params(0, 0, 0, 8/*xd->bd*/);
                        uint32_t list_index = (sequence_control_set_ptr->mrp_mode == 0) ? (best_rf_idx < 4 ? 0 : 1)
                                                                                        : (best_rf_idx < 2 ? 0 : 1);
                        uint32_t ref_pic_index = (sequence_control_set_ptr->mrp_mode == 0) ? (best_rf_idx >= 4 ? (best_rf_idx - 4) : best_rf_idx)
                                                                                           : (best_rf_idx >= 2 ? (best_rf_idx - 2) : best_rf_idx);
                        const int ref_mb_offset = mb_origin_y * input_picture_ptr->stride_y + mb_origin_x;
                        referenceObject = (EbReferenceObject*)picture_control_set_ptr->ref_pa_pic_ptr_array[list_index][ref_pic_index]->object_ptr;
                        ref_pic_ptr = /*is16bit ? (EbPictureBufferDesc*)referenceObject->reference_picture16bit : */(EbPictureBufferDesc*)referenceObject->reference_picture;
                        uint8_t *ref_mb = ref_pic_ptr->buffer_y + ref_mb_offset;
                        uint8_t *src_mb = input_picture_ptr->buffer_y + mb_origin_y * input_picture_ptr->stride_y + mb_origin_x;
                        MV best_mv = {x_curr_mv, y_curr_mv};
                        av1_build_inter_predictor(ref_mb, input_picture_ptr->stride_y, predictor, 16,
                                &best_mv, mb_origin_x, mb_origin_y, &inter_pred_params);
//printf("kelvin inter pred sb_index=%d, mb_origin_x=%d, mb_origin_y=%d\n", sb_index, mb_origin_x, mb_origin_y);
                    } else {
                        // intra recon
                        uint8_t *above_row;
                        uint8_t *left_col;
                        uint32_t mb_stride = (sequence_control_set_ptr->seq_header.max_frame_width + 15) / 16;
                        DECLARE_ALIGNED(16, uint8_t, left_data[MAX_TX_SIZE * 2 + 32]);
                        DECLARE_ALIGNED(16, uint8_t, above_data[MAX_TX_SIZE * 2 + 32]);

//printf("kelvin do intra pred sb_index=%d, mb_origin_x=%d, mb_origin_y=%d\n", sb_index, mb_origin_x, mb_origin_y);
                        above_row = above_data + 16;
                        left_col = left_data + 16;
                        TxSize tx_size = TX_16X16;
                        uint8_t *src = input_picture_ptr->buffer_y + picture_control_set_ptr->enhanced_picture_ptr->origin_x + mb_origin_x +
                            (picture_control_set_ptr->enhanced_picture_ptr->origin_y + mb_origin_y) * input_picture_ptr->stride_y;
#if USE_ORIGIN_YUV
                        update_neighbor_samples_array_open_loop(above_row - 1, left_col - 1, input_picture_ptr, input_picture_ptr->stride_y, mb_origin_x, mb_origin_y, 16, 16, picture_control_set_ptr);
#else
                        update_neighbor_samples_array_open_loop(above_row - 1, left_col - 1, input_picture_ptr, input_picture_ptr->stride_y, mb_origin_x, mb_origin_y, 16, 16);
#endif
                        uint8_t ois_intra_mode = ois_mb_results_ptr->intra_mode;
                        int32_t p_angle = av1_is_directional_mode((PredictionMode)ois_intra_mode) ? mode_to_angle_map[(PredictionMode)ois_intra_mode] : 0;
                        // Edge filter
                        if(av1_is_directional_mode((PredictionMode)ois_intra_mode) && 1/*sequence_control_set_ptr->seq_header.enable_intra_edge_filter*/) {
                            above_row = above_data + 16;
                            left_col  = left_data + 16;
                            filter_intra_edge(picture_control_set_ptr, ois_mb_results_ptr, ois_intra_mode, p_angle, mb_origin_x, mb_origin_y, above_row, left_col);
                        }
                        // PRED
                        intra_prediction_open_loop_mb(p_angle, ois_intra_mode, mb_origin_x, mb_origin_y, tx_size, above_row, left_col, predictor, 16/*dst_stride*/);
                    }

                    aom_subtract_block(16, 16, src_diff, 16, src_mb, input_picture_ptr->stride_y, predictor, 16);
                    wht_fwd_txfm(src_diff, 16, coeff, tx_size, 8/*xd->bd*/, 0/*s_cur_buf_hbd(xd)*/);

                    uint16_t eob;
                    get_quantize_error(&mb_plane, 0, coeff, qcoeff, dqcoeff, tx_size, &eob, &recon_error, &sse);

                    int rate_cost = rate_estimator(qcoeff, eob, tx_size);
                    if(eob) {
                        /*if(is16bit)
                            av1_inv_transform_recon16bit();
                        else*/
//printf("kelvincost1 poc%d sb_index=%d, mb_origin_xy=%d %d, before av1_inv_transform_recon8bit, eob=%d\n", picture_control_set_ptr->picture_number, sb_index, mb_origin_x, mb_origin_y, eob);
                            av1_inv_transform_recon8bit((int32_t*)dqcoeff, predictor, 16, predictor, 16, TX_16X16, DCT_DCT, PLANE_TYPE_Y, eob, 0 /*lossless*/);
                    }

                    ois_mb_results_ptr->recrf_dist = recon_error << (TPL_DEP_COST_SCALE_LOG2);
                    ois_mb_results_ptr->recrf_rate = rate_cost << TPL_DEP_COST_SCALE_LOG2;
                    if (best_mode != NEWMV) {
                        ois_mb_results_ptr->srcrf_dist = recon_error << (TPL_DEP_COST_SCALE_LOG2);
                        ois_mb_results_ptr->srcrf_rate = rate_cost << TPL_DEP_COST_SCALE_LOG2;
                    }
                    ois_mb_results_ptr->recrf_dist = AOMMAX(ois_mb_results_ptr->srcrf_dist, ois_mb_results_ptr->recrf_dist);
                    ois_mb_results_ptr->recrf_rate = AOMMAX(ois_mb_results_ptr->srcrf_rate, ois_mb_results_ptr->recrf_rate);

                    if (/*frame_idx && */best_rf_idx != -1) {
                        ois_mb_results_ptr->mv = final_best_mv;
                        ois_mb_results_ptr->ref_frame_poc = picture_control_set_ptr->ref_order_hint[best_rf_idx];
                    }

                    // Motion flow dependency dispenser.
                    result_model_store(picture_control_set_ptr, ois_mb_results_ptr, mb_origin_x, mb_origin_y, picture_width_in_mb);
                }
                pa_blk_index++;
            }
        }
    }

    return;
}

static int get_overlap_area(int grid_pos_row, int grid_pos_col, int ref_pos_row,
                            int ref_pos_col, int block, int/*BLOCK_SIZE*/ bsize) {
  int width = 0, height = 0;
  int bw = 4 << mi_size_wide_log2[bsize];
  int bh = 4 << mi_size_high_log2[bsize];

  switch (block) {
    case 0:
      width = grid_pos_col + bw - ref_pos_col;
      height = grid_pos_row + bh - ref_pos_row;
      break;
    case 1:
      width = ref_pos_col + bw - grid_pos_col;
      height = grid_pos_row + bh - ref_pos_row;
      break;
    case 2:
      width = grid_pos_col + bw - ref_pos_col;
      height = ref_pos_row + bh - grid_pos_row;
      break;
    case 3:
      width = ref_pos_col + bw - grid_pos_col;
      height = ref_pos_row + bh - grid_pos_row;
      break;
    default: assert(0);
  }

  return width * height;
}

static int round_floor(int ref_pos, int bsize_pix) {
  int round;
  if (ref_pos < 0)
    round = -(1 + (-ref_pos - 1) / bsize_pix);
  else
    round = ref_pos / bsize_pix;

  return round;
}

static int64_t delta_rate_cost(int64_t delta_rate, int64_t recrf_dist,
                               int64_t srcrf_dist, int pix_num) {
  double beta = (double)srcrf_dist / recrf_dist;
  int64_t rate_cost = delta_rate;

  if (srcrf_dist <= 128) return rate_cost;

  double dr =
      (double)(delta_rate >> (TPL_DEP_COST_SCALE_LOG2 + AV1_PROB_COST_SHIFT)) /
      pix_num;

  double log_den = log(beta) / log(2.0) + 2.0 * dr;

  if (log_den > log(10.0) / log(2.0)) {
    rate_cost = (int64_t)((log(1.0 / beta) * pix_num) / log(2.0) / 2.0);
    rate_cost <<= (TPL_DEP_COST_SCALE_LOG2 + AV1_PROB_COST_SHIFT);
    return rate_cost;
  }

  double num = pow(2.0, log_den);
  double den = num * beta + (1 - beta) * beta;

  rate_cost = (int64_t)((pix_num * log(num / den)) / log(2.0) / 2.0);

  rate_cost <<= (TPL_DEP_COST_SCALE_LOG2 + AV1_PROB_COST_SHIFT);

  return rate_cost;
}

static AOM_INLINE void tpl_model_update_b(PictureParentControlSet *ref_picture_control_set_ptr, OisMbResults *ois_mb_results_ptr,
                                          int mi_row, int mi_col,
                                          const int/*BLOCK_SIZE*/ bsize) {
  Av1Common *ref_cm = ref_picture_control_set_ptr->av1_cm;
  OisMbResults *ref_ois_mb_results_ptr;

  const int ref_pos_row = mi_row * MI_SIZE + (ois_mb_results_ptr->mv.row >> 3);
  const int ref_pos_col = mi_col * MI_SIZE + (ois_mb_results_ptr->mv.col >> 3);

  const int bw = 4 << mi_size_wide_log2[bsize];
  const int bh = 4 << mi_size_high_log2[bsize];
  const int mi_height = mi_size_high[bsize];
  const int mi_width = mi_size_wide[bsize];
  const int pix_num = bw * bh;

  // top-left on grid block location in pixel
  int grid_pos_row_base = round_floor(ref_pos_row, bh) * bh;
  int grid_pos_col_base = round_floor(ref_pos_col, bw) * bw;
  int block;

  int64_t cur_dep_dist = ois_mb_results_ptr->recrf_dist - ois_mb_results_ptr->srcrf_dist;
  int64_t mc_dep_dist = (int64_t)(
      ois_mb_results_ptr->mc_dep_dist *
      ((double)(ois_mb_results_ptr->recrf_dist - ois_mb_results_ptr->srcrf_dist) /
       ois_mb_results_ptr->recrf_dist));
  int64_t delta_rate = ois_mb_results_ptr->recrf_rate - ois_mb_results_ptr->srcrf_rate;
  int64_t mc_dep_rate =
      delta_rate_cost(ois_mb_results_ptr->mc_dep_rate, ois_mb_results_ptr->recrf_dist,
                      ois_mb_results_ptr->srcrf_dist, pix_num);

  for (block = 0; block < 4; ++block) {
    int grid_pos_row = grid_pos_row_base + bh * (block >> 1);
    int grid_pos_col = grid_pos_col_base + bw * (block & 0x01);

    if (grid_pos_row >= 0 && grid_pos_row < ref_cm->mi_rows * MI_SIZE &&
        grid_pos_col >= 0 && grid_pos_col < ref_cm->mi_cols * MI_SIZE) {
      int overlap_area = get_overlap_area(
          grid_pos_row, grid_pos_col, ref_pos_row, ref_pos_col, block, bsize);
      int ref_mi_row = round_floor(grid_pos_row, bh) * mi_height;
      int ref_mi_col = round_floor(grid_pos_col, bw) * mi_width;
      const int step = 1 << 2;//cpi->tpl_stats_block_mis_log2;

      for (int idy = 0; idy < mi_height; idy += step) {
        for (int idx = 0; idx < mi_width; idx += step) {
          ref_ois_mb_results_ptr = ref_picture_control_set_ptr->ois_mb_results[((ref_mi_row + idy) >> 2) * (ref_cm->mi_rows >> 2)  + ((ref_mi_col + idx) >> 2)]; //cpi->tpl_stats_block_mis_log2;
          ref_ois_mb_results_ptr->mc_dep_dist +=
              ((cur_dep_dist + mc_dep_dist) * overlap_area) / pix_num;
          ref_ois_mb_results_ptr->mc_dep_rate +=
              ((delta_rate + mc_dep_rate) * overlap_area) / pix_num;

          assert(overlap_area >= 0);
        }
      }
    }
  }
}

static AOM_INLINE void tpl_model_update(PictureParentControlSet *picture_control_set_array[60], int32_t frame_idx, int mi_row, int mi_col, const int/*BLOCK_SIZE*/ bsize, uint8_t frames_in_sw) {
  const int mi_height = mi_size_high[bsize];
  const int mi_width = mi_size_wide[bsize];
  const int step = 1 << 2; //cpi->tpl_stats_block_mis_log2;
  const int/*BLOCK_SIZE*/ block_size = BLOCK_16X16; //convert_length_to_bsize(MI_SIZE << cpi->tpl_stats_block_mis_log2);
  PictureParentControlSet *picture_control_set_ptr = picture_control_set_array[frame_idx];
  int i = 0;

  for (int idy = 0; idy < mi_height; idy += step) {
    for (int idx = 0; idx < mi_width; idx += step) {
      OisMbResults *ois_mb_results_ptr = picture_control_set_ptr->ois_mb_results[(((mi_row + idy) * mi_width) >> 4) + ((mi_col + idx) >> 2)];
      while(i<frames_in_sw && picture_control_set_array[i]->picture_number != ois_mb_results_ptr->ref_frame_poc)
        i++;
      //assert(picture_control_set_array[i]->picture_number < frame_idx); //kelvin to check
      if(i<frames_in_sw)
        tpl_model_update_b(picture_control_set_array[i], ois_mb_results_ptr, mi_row + idy, mi_col + idx, block_size);
    }
  }
}

void cutree_mc_flow_synthesizer(
    EncodeContext                   *encode_context_ptr,
    SequenceControlSet              *sequence_control_set_ptr,
    PictureParentControlSet         *picture_control_set_array[60],
    int32_t                          frame_idx,
    uint8_t                          frames_in_sw)
{
    Av1Common *cm = picture_control_set_array[frame_idx]->av1_cm;
    const int/*BLOCK_SIZE*/ bsize = BLOCK_16X16;
    const int mi_height = mi_size_high[bsize];
    const int mi_width = mi_size_wide[bsize];

    for (int mi_row = 0; mi_row < cm->mi_rows; mi_row += mi_height) {
        for (int mi_col = 0; mi_col < cm->mi_cols; mi_col += mi_width) {
            tpl_model_update(picture_control_set_array, frame_idx, mi_row, mi_col, bsize, frames_in_sw);
        }
    }
    return;
}

/************************************************
* Genrate CUTree MC Flow Synthesizer Based on Lookahead
** LAD Window: sliding window size
************************************************/
void update_mc_flow_synthesizer(
    EncodeContext                   *encode_context_ptr,
    SequenceControlSet              *sequence_control_set_ptr,
    PictureParentControlSet         *picture_control_set_ptr)
{
    InitialRateControlReorderEntry   *temporaryQueueEntryPtr;
    PictureParentControlSet          *temp_picture_control_set_ptr;
    PictureParentControlSet          *picture_control_set_array[60] = {NULL, };

    //uint32_t                         frame_idx_checked = 0;
    uint32_t                         inputQueueIndex;
    int32_t                          frame_idx, i;

    (void)sequence_control_set_ptr;

    // Walk the first N entries in the sliding window
    inputQueueIndex = encode_context_ptr->initial_rate_control_reorder_queue_head_index;
    for (frame_idx = 0; frame_idx < picture_control_set_ptr->frames_in_sw; frame_idx++) {
        temporaryQueueEntryPtr = encode_context_ptr->initial_rate_control_reorder_queue[inputQueueIndex];
        temp_picture_control_set_ptr = ((PictureParentControlSet*)(temporaryQueueEntryPtr->parent_pcs_wrapper_ptr)->object_ptr);

        //printf("kelvin ---> init_rc update_mc_flow_synthesizer curr picture_number=%d %d, inputQueueIndex=%d, temp picture_number=%d, temp decode_order=%d, frames_in_sw=%d\n", picture_control_set_ptr->picture_number, frame_idx, inputQueueIndex, temp_picture_control_set_ptr->picture_number, temp_picture_control_set_ptr->decode_order, picture_control_set_ptr->frames_in_sw);

        // sort to be decode order
        if(frame_idx == 0) {
            picture_control_set_array[0] = temp_picture_control_set_ptr;
        } else {
            for (i = 0; i < frame_idx; i++) {
                if (temp_picture_control_set_ptr->decode_order < picture_control_set_array[i]->decode_order) {
                    for (int32_t j = frame_idx; j > i; j--)
                        picture_control_set_array[j] = picture_control_set_array[j-1];
                    picture_control_set_array[i] = temp_picture_control_set_ptr;
                    break;
                }
            }
            if (i == frame_idx)
                picture_control_set_array[i] = temp_picture_control_set_ptr;
        }

        // Increment the inputQueueIndex Iterator
        inputQueueIndex = (inputQueueIndex == INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1) ? 0 : inputQueueIndex + 1;
    }

    for(frame_idx = picture_control_set_ptr->frames_in_sw - 1; frame_idx > 0; frame_idx--) {
        printf("kelvin ---> init_rc update_mc_flow_synthesizer frame_idx=%d, reordered picture_number=%d, decode_order=%d\n", frame_idx, picture_control_set_array[frame_idx]->picture_number, picture_control_set_array[frame_idx]->decode_order);
        //kelvinhack
        cutree_mc_flow_synthesizer(encode_context_ptr, sequence_control_set_ptr, picture_control_set_array, frame_idx, picture_control_set_ptr->frames_in_sw); // in decode order
    }

    return;
}

#endif
/* Initial Rate Control Kernel */

/*********************************************************************************
*
* @brief
*  The Initial Rate Control process determines the initial bit budget for each picture
*  depending on the data gathered in the Picture Analysis and Motion Estimation processes
*  as well as the settings determined in the Picture Decision process.
*
* @par Description:
*  The Initial Rate Control process employs a sliding window buffer to analyze
*  multiple pictures if a delay is allowed. Note that no reference picture data is
*  used in this process.
*
* @param[in] Picture
*  The Initial Rate Control Kernel takes a picture and determines the initial bit budget
*  for each picture depending on the data that was gathered in Picture Analysis and
*  Motion Estimation processes
*
* @param[out] Bit Budget
*  Bit Budget is the amount of budgetted bits for a picture
*
* @remarks
*  Temporal noise reduction is currently performed in Initial Rate Control Process.
*  In the future we might decide to move it to Motion Analysis Process.
*
********************************************************************************/
void *initial_rate_control_kernel(void *input_ptr) {
    EbThreadContext *          thread_context_ptr = (EbThreadContext *)input_ptr;
    InitialRateControlContext *context_ptr = (InitialRateControlContext *)thread_context_ptr->priv;
    PictureParentControlSet *  pcs_ptr;
    PictureParentControlSet *  pcs_ptr_temp;
    EncodeContext *            encode_context_ptr;
    SequenceControlSet *       scs_ptr;

    EbObjectWrapper *        in_results_wrapper_ptr;
    MotionEstimationResults *in_results_ptr;

    EbObjectWrapper *          out_results_wrapper_ptr;
    InitialRateControlResults *out_results_ptr;

    // Queue variables
    uint32_t                        queue_entry_index_temp;
    uint32_t                        queue_entry_index_temp2;
    InitialRateControlReorderEntry *queue_entry_ptr;

    EbBool           move_slide_window_flag = EB_TRUE;
    EbBool           end_of_sequence_flag   = EB_TRUE;
    uint8_t          frames_in_sw;
    uint8_t          temporal_layer_index;
    EbObjectWrapper *reference_picture_wrapper_ptr;

    // Segments
    uint32_t segment_index;
    for (;;) {
        // Get Input Full Object
        eb_get_full_object(context_ptr->motion_estimation_results_input_fifo_ptr,
                           &in_results_wrapper_ptr);

        in_results_ptr = (MotionEstimationResults *)in_results_wrapper_ptr->object_ptr;
        pcs_ptr        = (PictureParentControlSet *)in_results_ptr->pcs_wrapper_ptr->object_ptr;

        segment_index = in_results_ptr->segment_index;

        // Set the segment mask
        SEGMENT_COMPLETION_MASK_SET(pcs_ptr->me_segments_completion_mask, segment_index);

        // If the picture is complete, proceed
        if (SEGMENT_COMPLETION_MASK_TEST(pcs_ptr->me_segments_completion_mask,
                                         pcs_ptr->me_segments_total_count)) {
            scs_ptr            = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
            encode_context_ptr = (EncodeContext *)scs_ptr->encode_context_ptr;
            // Mark picture when global motion is detected using ME results
            //reset intra_coded_estimation_sb
            me_based_global_motion_detection(pcs_ptr);
#if CUTREE_LA
            if(picture_control_set_ptr->picture_number == 0)
                printf("kelvin ---> input sqc look_ahead_distance=%d, enable_cutree_in_la=%d\n", sequence_control_set_ptr->static_config.look_ahead_distance, sequence_control_set_ptr->static_config.enable_cutree_in_la);
            if (sequence_control_set_ptr->static_config.look_ahead_distance == 0 || sequence_control_set_ptr->static_config.enable_cutree_in_la == 0) {
                // Release Pa Ref pictures when not needed
                ReleasePaReferenceObjects(
                    sequence_control_set_ptr,
                    picture_control_set_ptr);
            }
            if (sequence_control_set_ptr->static_config.look_ahead_distance != 0 && sequence_control_set_ptr->static_config.enable_cutree_in_la) {
                //if (picture_control_set_ptr->slice_type != I_SLICE)
                    //kelvinhack
                    cutree_mc_flow_dispenser(encode_context_ptr, sequence_control_set_ptr, picture_control_set_ptr);
            }
            printf("kelvin ---> init_rc picture_number=%d to releasePaReference, look_ahead_distance=%d\n", picture_control_set_ptr->picture_number, sequence_control_set_ptr->static_config.look_ahead_distance);
#else
            // Release Pa Ref pictures when not needed
            release_pa_reference_objects(scs_ptr, pcs_ptr);
#endif

            //****************************************************
            // Picture resizing for super-res tool
            //****************************************************

            // Scale picture if super-res is used
            if(scs_ptr->static_config.superres_mode > SUPERRES_NONE){
                init_resize_picture(pcs_ptr->scs_ptr,
                                    pcs_ptr);
            }

            //****************************************************
            // Input Motion Analysis Results into Reordering Queue
            //****************************************************

            if (!pcs_ptr->is_overlay)
                // Determine offset from the Head Ptr
                queue_entry_ptr =
                    determine_picture_offset_in_queue(encode_context_ptr, pcs_ptr, in_results_ptr);

            if (scs_ptr->static_config.rate_control_mode) {
                if (scs_ptr->static_config.look_ahead_distance != 0) {
                    // Getting the Histogram Queue Data
                    get_histogram_queue_data(scs_ptr, encode_context_ptr, pcs_ptr);
                }
            }

            for (temporal_layer_index = 0; temporal_layer_index < EB_MAX_TEMPORAL_LAYERS;
                 temporal_layer_index++)
                pcs_ptr->frames_in_interval[temporal_layer_index] = 0;
            pcs_ptr->frames_in_sw          = 0;
            pcs_ptr->historgram_life_count = 0;
            pcs_ptr->scene_change_in_gop   = EB_FALSE;
            move_slide_window_flag = EB_TRUE;
            while (move_slide_window_flag) {
                // Check if the sliding window condition is valid
                queue_entry_index_temp =
                    encode_context_ptr->initial_rate_control_reorder_queue_head_index;
                if (encode_context_ptr->initial_rate_control_reorder_queue[queue_entry_index_temp]
                        ->parent_pcs_wrapper_ptr != EB_NULL)
                    end_of_sequence_flag =
                        (((PictureParentControlSet
                               *)(encode_context_ptr
                                      ->initial_rate_control_reorder_queue[queue_entry_index_temp]
                                      ->parent_pcs_wrapper_ptr)
                              ->object_ptr))
                            ->end_of_sequence_flag;
                else
                    end_of_sequence_flag = EB_FALSE;
                frames_in_sw = 0;
                while (move_slide_window_flag && !end_of_sequence_flag &&
                       queue_entry_index_temp <=
                           encode_context_ptr->initial_rate_control_reorder_queue_head_index +
                               scs_ptr->static_config.look_ahead_distance) {
                    // frames_in_sw <= scs_ptr->static_config.look_ahead_distance){
                    frames_in_sw++;

                    queue_entry_index_temp2 =
                        (queue_entry_index_temp > INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                            ? queue_entry_index_temp - INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH
                            : queue_entry_index_temp;

                    move_slide_window_flag =
                        (EbBool)(move_slide_window_flag &&
                                 (encode_context_ptr
                                      ->initial_rate_control_reorder_queue[queue_entry_index_temp2]
                                      ->parent_pcs_wrapper_ptr != EB_NULL));
                    if (encode_context_ptr
                            ->initial_rate_control_reorder_queue[queue_entry_index_temp2]
                            ->parent_pcs_wrapper_ptr != EB_NULL) {
                        // check if it is the last frame. If we have reached the last frame, we would output the buffered frames in the Queue.
                        end_of_sequence_flag =
                            ((PictureParentControlSet *)(encode_context_ptr
                                                             ->initial_rate_control_reorder_queue
                                                                 [queue_entry_index_temp2]
                                                             ->parent_pcs_wrapper_ptr)
                                 ->object_ptr)
                                ->end_of_sequence_flag;
                    } else
                        end_of_sequence_flag = EB_FALSE;
                    queue_entry_index_temp++;
                }

                if (move_slide_window_flag) {
                    //get a new entry spot
                    queue_entry_ptr =
                        encode_context_ptr->initial_rate_control_reorder_queue
                            [encode_context_ptr->initial_rate_control_reorder_queue_head_index];
                    pcs_ptr = ((PictureParentControlSet *)(queue_entry_ptr->parent_pcs_wrapper_ptr)
                                   ->object_ptr);
                    scs_ptr = (SequenceControlSet *)pcs_ptr->scs_wrapper_ptr->object_ptr;
                    // overlay picture was not added to the queue. For the alt_ref picture with an overlay picture, it loops on both alt ref and overlay pictures
                    uint8_t has_overlay = pcs_ptr->is_alt_ref ? 1 : 0;
                    for (uint8_t loop_index = 0; loop_index <= has_overlay; loop_index++) {
                        if (loop_index) pcs_ptr = pcs_ptr->overlay_ppcs_ptr;
                        pcs_ptr->frames_in_sw = frames_in_sw;
                        queue_entry_index_temp =
                            encode_context_ptr->initial_rate_control_reorder_queue_head_index;
                        end_of_sequence_flag = EB_FALSE;
                        // find the frames_in_interval for the peroid I frames
                        while (
                            !end_of_sequence_flag &&
                            queue_entry_index_temp <=
                                encode_context_ptr->initial_rate_control_reorder_queue_head_index +
                                    scs_ptr->static_config.look_ahead_distance) {
                            queue_entry_index_temp2 =
                                (queue_entry_index_temp >
                                 INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                                    ? queue_entry_index_temp -
                                          INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH
                                    : queue_entry_index_temp;
                            pcs_ptr_temp = ((PictureParentControlSet
                                                 *)(encode_context_ptr
                                                        ->initial_rate_control_reorder_queue
                                                            [queue_entry_index_temp2]
                                                        ->parent_pcs_wrapper_ptr)
                                                ->object_ptr);
                            if (scs_ptr->intra_period_length != -1) {
                                if (pcs_ptr->picture_number %
                                        ((scs_ptr->intra_period_length + 1)) ==
                                    0) {
                                    pcs_ptr
                                        ->frames_in_interval[pcs_ptr_temp->temporal_layer_index]++;
                                    if (pcs_ptr_temp->scene_change_flag)
                                        pcs_ptr->scene_change_in_gop = EB_TRUE;
                                }
                            }

                            pcs_ptr_temp->historgram_life_count++;
                            end_of_sequence_flag = pcs_ptr_temp->end_of_sequence_flag;
                            queue_entry_index_temp++;
                        }

                        if ((scs_ptr->static_config.look_ahead_distance != 0) &&
                            (frames_in_sw < (scs_ptr->static_config.look_ahead_distance + 1)))
                            pcs_ptr->end_of_sequence_region = EB_TRUE;
                        else
                            pcs_ptr->end_of_sequence_region = EB_FALSE;

                        if (scs_ptr->static_config.rate_control_mode) {
                            // Determine offset from the Head Ptr for HLRC histogram queue and set the life count
                            if (scs_ptr->static_config.look_ahead_distance != 0) {
                                // Update Histogram Queue Entry Life count
                                update_histogram_queue_entry(
                                    scs_ptr, encode_context_ptr, pcs_ptr, frames_in_sw);
                            }
                        }

                        // Mark each input picture as PAN or not
                        // If a lookahead is present then check PAN for a period of time
                        if (!pcs_ptr->end_of_sequence_flag &&
                            scs_ptr->static_config.look_ahead_distance != 0) {
                            // Check for Pan,Tilt, Zoom and other global motion detectors over the future pictures in the lookahead
                            update_global_motion_detection_over_time(
                                encode_context_ptr, scs_ptr, pcs_ptr);
                        } else {
                            if (pcs_ptr->slice_type != I_SLICE) detect_global_motion(pcs_ptr);
                        }

                        // BACKGROUND ENHANCEMENT Part II
                        if (!pcs_ptr->end_of_sequence_flag &&
                            scs_ptr->static_config.look_ahead_distance != 0) {
                            // Update BEA information based on Lookahead information
                            update_bea_info_over_time(encode_context_ptr, pcs_ptr);
                        } else {
                            // Reset zzCost information to default When there's no lookahead available
                            init_zz_cost_info(pcs_ptr);
                        }

                        // Use the temporal layer 0 is_sb_motion_field_non_uniform array for all the other layer pictures in the mini GOP
                        if (!pcs_ptr->end_of_sequence_flag &&
                            scs_ptr->static_config.look_ahead_distance != 0) {
                            // Updat uniformly moving SBs based on Collocated SBs in LookAhead window
                            update_motion_field_uniformity_over_time(
                                encode_context_ptr, scs_ptr, pcs_ptr);
                        }
                        // Get Empty Reference Picture Object
                        eb_get_empty_object(
                            scs_ptr->encode_context_ptr->reference_picture_pool_fifo_ptr,
                            &reference_picture_wrapper_ptr);
                        if (loop_index) {
                            pcs_ptr->reference_picture_wrapper_ptr = reference_picture_wrapper_ptr;
                            // Give the new Reference a nominal live_count of 1
                            eb_object_inc_live_count(pcs_ptr->reference_picture_wrapper_ptr, 1);
                        } else {
                            ((PictureParentControlSet *)(queue_entry_ptr->parent_pcs_wrapper_ptr
                                                             ->object_ptr))
                                ->reference_picture_wrapper_ptr = reference_picture_wrapper_ptr;
                            // Give the new Reference a nominal live_count of 1
                            eb_object_inc_live_count(
                                ((PictureParentControlSet *)(queue_entry_ptr->parent_pcs_wrapper_ptr
                                                                 ->object_ptr))
                                    ->reference_picture_wrapper_ptr,
                                1);
                        }
                        pcs_ptr->stat_struct_first_pass_ptr =
                            pcs_ptr->is_used_as_reference_flag
                                ? &((EbReferenceObject *)
                                        pcs_ptr->reference_picture_wrapper_ptr->object_ptr)
                                       ->stat_struct
                                : &pcs_ptr->stat_struct;
                        if (scs_ptr->use_output_stat_file)
                            memset(pcs_ptr->stat_struct_first_pass_ptr, 0, sizeof(StatStruct));
#if CUTREE_LA
                        if (sequence_control_set_ptr->static_config.look_ahead_distance != 0 &&
                            sequence_control_set_ptr->static_config.enable_cutree_in_la &&
                            picture_control_set_ptr->temporal_layer_index == 0) {
                            update_mc_flow_synthesizer(encode_context_ptr, sequence_control_set_ptr, picture_control_set_ptr);
                        }
#endif

                        // Get Empty Results Object
                        eb_get_empty_object(
                            context_ptr->initialrate_control_results_output_fifo_ptr,
                            &out_results_wrapper_ptr);

                        out_results_ptr =
                            (InitialRateControlResults *)out_results_wrapper_ptr->object_ptr;

                        if (loop_index)
                            out_results_ptr->pcs_wrapper_ptr = pcs_ptr->p_pcs_wrapper_ptr;
                        else
                            out_results_ptr->pcs_wrapper_ptr =
                                queue_entry_ptr->parent_pcs_wrapper_ptr;
#if CUTREE_LA
printf("kelvin ---> loop_index=%d, output picture_number=%d\n", loop_index, loop_index ? picture_control_set_ptr->picture_number : queueEntryPtr->picture_number );
                        if (sequence_control_set_ptr->static_config.look_ahead_distance != 0 && sequence_control_set_ptr->static_config.enable_cutree_in_la
                            && ((has_overlay == 0 && loop_index == 0) || (has_overlay == 1 && loop_index == 1))) {
                            // Release Pa Ref pictures when not needed
                            ReleasePaReferenceObjects(
                                sequence_control_set_ptr,
                                picture_control_set_ptr);
                                //loop_index ? picture_control_set_ptr : queueEntryPtr);
                        }
#endif
                        // Post the Full Results Object
                        eb_post_full_object(out_results_wrapper_ptr);
                    }
                    // Reset the Reorder Queue Entry
                    queue_entry_ptr->picture_number += INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH;
                    queue_entry_ptr->parent_pcs_wrapper_ptr = (EbObjectWrapper *)EB_NULL;

                    // Increment the Reorder Queue head Ptr
                    encode_context_ptr->initial_rate_control_reorder_queue_head_index =
                        (encode_context_ptr->initial_rate_control_reorder_queue_head_index ==
                         INITIAL_RATE_CONTROL_REORDER_QUEUE_MAX_DEPTH - 1)
                            ? 0
                            : encode_context_ptr->initial_rate_control_reorder_queue_head_index + 1;

                    queue_entry_ptr =
                        encode_context_ptr->initial_rate_control_reorder_queue
                            [encode_context_ptr->initial_rate_control_reorder_queue_head_index];
                }
            }
        }

        // Release the Input Results
        eb_release_object(in_results_wrapper_ptr);
    }
    return EB_NULL;
}
