#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <assert.h>

#include "nalu.h"

// 从SPS中解析宽高
void h264_parse_sps_wh(const sps_t *sps, int *width, int *height);
// 从SPS中解析帧率，可能返回0，表示该SPS未携带帧率信息
int h264_parse_sps_fps(const sps_t *sps);

void h264_parse_sps_wh(const sps_t *sps, int *width, int *height) {
    *width = (sps->pic_width_in_mbs_minus1 + 1) * 16;
    *height = (sps->pic_height_in_map_units_minus1 + 1) * 16;

    if(sps->frame_cropping_flag) {
        unsigned int crop_unit_x;
        unsigned int crop_unit_y;
        switch (sps->chroma_format_idc)
        {
        case YUV_4_0_0:
            crop_unit_x = 1;
            crop_unit_y = 2 - sps->frame_mbs_only_flag;
            break;
        case YUV_4_2_0:
            crop_unit_x = 2;
            crop_unit_y = 2 * (2 - sps->frame_mbs_only_flag);
            break;
        case YUV_4_2_2:
            crop_unit_x = 2;
            crop_unit_y = 2 - sps->frame_mbs_only_flag;
            break;
        case YUV_4_4_4:
            crop_unit_x = 1;
            crop_unit_y = 2 - sps->frame_mbs_only_flag;
            break;
        default:
            printf("Invalid sps!!!\n");
            return;
        }
        *width -= crop_unit_x * (sps->frame_crop_left_offset + sps->frame_crop_right_offset);
        *height -= crop_unit_y * (sps->frame_crop_top_offset + sps->frame_crop_bottom_offset);
    }
}

int h264_parse_sps_fps(const sps_t *sps) {
    if(sps->vui_parameters_present_flag && sps->vui_parameters.timing_info_present_flag) {
        int fps = sps->vui_parameters.time_scale / sps->vui_parameters.num_units_in_tick;
        fps /= 2;
        return fps;
    } else {
        return 0;
    }
}

int main(int argc, char *argv[]) {
    if(argc < 2) {
        printf("Usage: %s xx.h264\n", argv[0]);
        return -1;
    }

    FILE *fp = NULL;
    int filesize = 0;
    uint8_t *buf = NULL;

    fp = fopen(argv[1], "rb");
    assert(fp);
    
    fseek(fp, 0, SEEK_END);
    filesize = (int)ftell(fp);
    printf("filesize: %d\n", filesize);

    buf = (uint8_t*)malloc(filesize);
    assert(buf);

    rewind(fp);
    fread(buf, sizeof(uint8_t), filesize, fp);
    fclose(fp);

    int nalu_i = 0;
    int nalu_size = 0;
    int cur_nal_start = -1;
    int cur_nal_end = -1;

    while((nalu_size = find_nal_unit(buf, filesize, &cur_nal_start, &cur_nal_end)) > 0) {
        printf("nalu: %03d, start: 0x%08x, end: 0x%08x, size: %-6d", nalu_i++, cur_nal_start, cur_nal_end, nalu_size);
        
        nal_t *nal = allocNAL();
        assert(nal);

        read_and_parse_nal_unit(nal, &buf[cur_nal_start], nalu_size);
        if(nal->nal_unit_type == NALU_TYPE_SPS) {
            int width, height;
            h264_parse_sps_wh(&nal->sps, &width, &height);
            printf("width: %d, height: %d, ", width, height);
            printf("fps: %d", h264_parse_sps_fps(&nal->sps));
        }
        printf("\n");
        
        freeNAL(nal);
    }

    free(buf);

    return 0;
}