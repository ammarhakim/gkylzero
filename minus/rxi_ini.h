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
 * - Adding API documentation
 */

#pragma once

typedef struct rxi_ini_t rxi_ini_t;

/**
 * Read file "filename" and return ini reader object. NULL if read
 * failed.
 *
 * @param filename Name of file to return
 * @return pointer to ini reader. NULL if read failed
 */
rxi_ini_t* rxi_ini_load(const char *filename);

/**
 * Free the reader when done with it.
 *
 * @param ini Reader to free
 */
void rxi_ini_free(rxi_ini_t *ini);

/**
 * Value associated with "key" in "section". Returns NULL if such
 * section/key does not exist.
 *
 * @param ini Ini reader
 * @param section Name of section to read from
 * @param key Key to read
 * @return Value associated with key.
 */
const char* rxi_ini_get(rxi_ini_t *ini, const char *section, const char *key);

/**
 * Read value into provided output variables specified by "scanfmt"
 * with "key" in "section". Returns 1 if read succeeded and 0
 * otherwise. The scanfmt is the one used by formatted output family
 * of function in C. Please see table at
 * https://en.cppreference.com/w/c/io/fprintf for full format
 * specifications.
 *
 * NOTE: You must use the *exact* format specifier or you will get
 * junk. For example, to read a double "charge" from [electron] and
 * store it in elc_charge do:
 *
 * rxi_ini_sget(ini, "electron", "charge", "%lg", &elc_charge);
 *
 * The format-specifier "%g" will NOT work to read double as it
 * specifies float and not double.
 *
 * @param ini Ini reader
 * @param section Name of section to read from
 * @param key Key to read
 * @param scanfmt Format specified strong. See example above
 * @return 0 if read failed and 1 otherwise
 */
int rxi_ini_sget(rxi_ini_t *ini, const char *section, const char *key, const char *scanfmt, ...);
