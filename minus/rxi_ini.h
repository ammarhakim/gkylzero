/**
 * Copyright (c) 2016 rxi
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the MIT license. See `ini.c` for details.
 * 
 * Modifications by Ammar Hakim, April 2021:
 *
 * - Appending rxi to public API
 * - Changed ini_sget to take variable number of arguments
 */

#pragma once

#define INI_VERSION "0.1.1"

typedef struct rxi_ini_t rxi_ini_t;

rxi_ini_t* rxi_ini_load(const char *filename);
void rxi_ini_free(rxi_ini_t *ini);
int rxi_num_sections(rxi_ini_t *ini);
const char* rxi_ini_get(rxi_ini_t *ini, const char *section, const char *key);
int rxi_ini_sget(rxi_ini_t *ini, const char *section, const char *key, const char *scanfmt, ...);
