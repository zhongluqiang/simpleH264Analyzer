#ifndef __NALU_H__
#define __NALU_H__

#include <stdint.h>
#include "typedef.h"

nal_t *allocNAL(void);
void freeNAL(nal_t *nal);

int find_nal_unit(const uint8_t *buf, int size, int *cur_nal_start, int *cur_nal_end);
int read_and_parse_nal_unit(nal_t *nal, uint8_t *buf, int size);

#endif