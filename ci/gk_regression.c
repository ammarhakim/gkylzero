#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

int system(const char *command);

int
main(int argc, char **argv)
{
  system("clear");
  system("mkdir -p output");

  printf("** Gkeyll Gyrokinetics Automated Regression System **\n");
  printf("Please select an option to proceed:\n\n");
  printf("1 - Run Full Regression Suite\n");
  printf("2 - View Regression Results\n");
  printf("3 - (Re)generate Accepted Results\n");
  printf("4 - Run Specific Test\n");
  printf("5 - (Re)generate Specific Accepted Result\n");

  int option;
  scanf("%d", &option);
  printf("\n");

  if (option == 1) {
    int sheath_1x2v_counter = 0;

    FILE *sheath_1x2v_counter_ptr = fopen("output/sheath_1x2v_counter.dat", "r");
    if (sheath_1x2v_counter_ptr != NULL) {
      fscanf(sheath_1x2v_counter_ptr, "%d", &sheath_1x2v_counter);
    }
    fclose(sheath_1x2v_counter_ptr);

    sheath_1x2v_counter += 1;

    sheath_1x2v_counter_ptr = fopen("output/sheath_1x2v_counter.dat", "w");
    fprintf(sheath_1x2v_counter_ptr, "%d", sheath_1x2v_counter);
    fclose(sheath_1x2v_counter_ptr);

    printf("Running 1x2v Sheath Boundary Test with p = 1...\n");
    system("cd ../; rm -rf ./gk_sheath_1x2v_p1-stat.json");
    system("cd ../; make build/regression/rt_gk_sheath_1x2v_p1 > /dev/null 2>&1");
    char sheath_1x2v_buffer[256];
    snprintf(sheath_1x2v_buffer, 256, "cd ../; ./build/regression/rt_gk_sheath_1x2v_p1 -m > ./ci/output/rt_gk_sheath_1x2v_p1_%d.dat 2>&1", sheath_1x2v_counter);
    system(sheath_1x2v_buffer);
    char sheath_1x2v_buffer2[256];
    snprintf(sheath_1x2v_buffer2, 256, "cd ../; mv ./gk_sheath_1x2v_p1-stat.json ci/output/gk_sheath_1x2v_p1-stat_%d.json", sheath_1x2v_counter);
    system(sheath_1x2v_buffer2);
    printf("Finished 1x2v Sheath Boundary Test with p = 1.\n\n");

    int sheath_2x2v_counter = 0;

    FILE *sheath_2x2v_counter_ptr = fopen("output/sheath_2x2v_counter.dat", "r");
    if (sheath_2x2v_counter_ptr != NULL) {
      fscanf(sheath_2x2v_counter_ptr, "%d", &sheath_2x2v_counter);
    }
    fclose(sheath_2x2v_counter_ptr);

    sheath_2x2v_counter += 1;

    sheath_2x2v_counter_ptr = fopen("output/sheath_2x2v_counter.dat", "w");
    fprintf(sheath_2x2v_counter_ptr, "%d", sheath_2x2v_counter);
    fclose(sheath_2x2v_counter_ptr);

    printf("Running 2x2v Sheath Boundary Test with p = 1...\n");
    system("cd ../; rm -rf ./gk_sheath_2x2v_p1-stat.json");
    system("cd ../; make build/regression/rt_gk_sheath_2x2v_p1 > /dev/null 2>&1");
    char sheath_2x2v_buffer[256];
    snprintf(sheath_2x2v_buffer, 256, "cd ../; ./build/regression/rt_gk_sheath_2x2v_p1 -m > ./ci/output/rt_gk_sheath_2x2v_p1_%d.dat 2>&1", sheath_2x2v_counter);
    system(sheath_2x2v_buffer);
    char sheath_2x2v_buffer2[256];
    snprintf(sheath_2x2v_buffer2, 256, "cd ../; mv ./gk_sheath_2x2v_p1-stat.json ci/output/gk_sheath_2x2v_p1-stat_%d.json", sheath_2x2v_counter);
    system(sheath_2x2v_buffer2);
    printf("Finished 2x2v Sheath Boundary Test with p = 1.\n\n");

    int sheath_3x2v_counter = 0;

    FILE *sheath_3x2v_counter_ptr = fopen("output/sheath_3x2v_counter.dat", "r");
    if (sheath_3x2v_counter_ptr != NULL) {
      fscanf(sheath_3x2v_counter_ptr, "%d", &sheath_3x2v_counter);
    }
    fclose(sheath_3x2v_counter_ptr);

    sheath_3x2v_counter += 1;

    sheath_3x2v_counter_ptr = fopen("output/sheath_3x2v_counter.dat", "w");
    fprintf(sheath_3x2v_counter_ptr, "%d", sheath_3x2v_counter);
    fclose(sheath_3x2v_counter_ptr);

    printf("Running 3x2v Sheath Boundary Test with p = 1...\n");
    system("cd ../; rm -rf ./gk_sheath_3x2v_p1-stat.json");
    system("cd ../; make build/regression/rt_gk_sheath_3x2v_p1 > /dev/null 2>&1");
    char sheath_3x2v_buffer[256];
    snprintf(sheath_3x2v_buffer, 256, "cd ../; ./build/regression/rt_gk_sheath_3x2v_p1 -m > ./ci/output/rt_gk_sheath_3x2v_p1_%d.dat 2>&1", sheath_3x2v_counter);
    system(sheath_3x2v_buffer);
    char sheath_3x2v_buffer2[256];
    snprintf(sheath_3x2v_buffer2, 256, "cd ../; mv ./gk_sheath_3x2v_p1-stat.json ci/output/gk_sheath_3x2v_p1-stat_%d.json", sheath_3x2v_counter);
    system(sheath_3x2v_buffer2);
    printf("Finished 3x2v Sheath Boundary Test with p = 1.\n\n");

    int lapd_cart_3x2v_counter = 0;

    FILE *lapd_cart_3x2v_counter_ptr = fopen("output/lapd_cart_3x2v_counter.dat", "r");
    if (lapd_cart_3x2v_counter_ptr != NULL) {
      fscanf(lapd_cart_3x2v_counter_ptr, "%d", &lapd_cart_3x2v_counter);
    }
    fclose(lapd_cart_3x2v_counter_ptr);

    lapd_cart_3x2v_counter += 1;

    lapd_cart_3x2v_counter_ptr = fopen("output/lapd_cart_3x2v_counter.dat", "w");
    fprintf(lapd_cart_3x2v_counter_ptr, "%d", lapd_cart_3x2v_counter);
    fclose(lapd_cart_3x2v_counter_ptr);

    printf("Running 3x2v LAPD Test (in Cartesian coordinates) with p = 1...\n");
    system("cd ../; rm -rf ./gk_lapd_cart_3x2v_p1-stat.json");
    system("cd ../; make build/regression/rt_gk_lapd_cart_3x2v_p1 > /dev/null 2>&1");
    char lapd_cart_3x2v_buffer[256];
    snprintf(lapd_cart_3x2v_buffer, 256, "cd ../; ./build/regression/rt_gk_lapd_cart_3x2v_p1 -m > ./ci/output/rt_gk_lapd_cart_3x2v_p1_%d.dat 2>&1", lapd_cart_3x2v_counter);
    system(lapd_cart_3x2v_buffer);
    char lapd_cart_3x2v_buffer2[256];
    snprintf(lapd_cart_3x2v_buffer2, 256, "cd ../; mv ./gk_lapd_cart_3x2v_p1-stat.json ci/output/gk_lapd_cart_3x2v_p1-stat_%d.json", lapd_cart_3x2v_counter);
    system(lapd_cart_3x2v_buffer2);
    printf("Finished 3x2v LAPD Test (in Cartesian coordinates) with p = 1.\n\n");

    int lapd_cyl_3x2v_counter = 0;

    FILE *lapd_cyl_3x2v_counter_ptr = fopen("output/lapd_cyl_3x2v_counter.dat", "r");
    if (lapd_cyl_3x2v_counter_ptr != NULL) {
      fscanf(lapd_cyl_3x2v_counter_ptr, "%d", &lapd_cyl_3x2v_counter);
    }
    fclose(lapd_cyl_3x2v_counter_ptr);

    lapd_cyl_3x2v_counter += 1;

    lapd_cyl_3x2v_counter_ptr = fopen("output/lapd_cyl_3x2v_counter.dat", "w");
    fprintf(lapd_cyl_3x2v_counter_ptr, "%d", lapd_cyl_3x2v_counter);
    fclose(lapd_cyl_3x2v_counter_ptr);

    printf("Running 3x2v LAPD Test (in cylindrical coordinates) with p = 1...\n");
    system("cd ../; rm -rf ./gk_lapd_cyl_3x2v_p1-stat.json");
    system("cd ../; make build/regression/rt_gk_lapd_cyl_3x2v_p1 > /dev/null 2>&1");
    char lapd_cyl_3x2v_buffer[256];
    snprintf(lapd_cyl_3x2v_buffer, 256, "cd ../; ./build/regression/rt_gk_lapd_cyl_3x2v_p1 -m > ./ci/output/rt_gk_lapd_cyl_3x2v_p1_%d.dat 2>&1", lapd_cyl_3x2v_counter);
    system(lapd_cyl_3x2v_buffer);
    char lapd_cyl_3x2v_buffer2[256];
    snprintf(lapd_cyl_3x2v_buffer2, 256, "cd ../; mv ./gk_lapd_cyl_3x2v_p1-stat.json ci/output/gk_lapd_cyl_3x2v_p1-stat_%d.json", lapd_cyl_3x2v_counter);
    system(lapd_cyl_3x2v_buffer2);
    printf("Finished 3x2v LAPD Test (in cylindrical coordinates) with p = 1.\n\n");
  }
  else if (option == 2) {
    printf("Please select the test whose results you wish to view:\n\n");
    printf("1 - 1x2v Sheath Boundary Test with p = 1\n");
    printf("2 - 2x2v Sheath Boundary Test with p = 1\n");
    printf("3 - 3x3v Sheath Boundary Test with p = 1\n");
    printf("4 - 3x2v LAPD Test (in Cartesian coordinates) with p = 1\n");
    printf("5 - 3x2v LAPD Test (in cylindrical coordinates) with p = 1\n");

    int option2;
    scanf("%d", &option2);
    printf("\n");

    if (option2 == 1) {
      printf("1x2v Sheath Boundary Test with p = 1:\n\n");

      int sheath_1x2v_counter = 0;

      FILE *sheath_1x2v_counter_ptr = fopen("output/sheath_1x2v_counter.dat", "r");
      if (sheath_1x2v_counter_ptr != NULL) {
        fscanf(sheath_1x2v_counter_ptr, "%d", &sheath_1x2v_counter);
      }
      fclose(sheath_1x2v_counter_ptr);

      int sheath_1x2v_updatecalls[sheath_1x2v_counter + 1];
      int sheath_1x2v_forwardeuler[sheath_1x2v_counter + 1];
      int sheath_1x2v_rk2failures[sheath_1x2v_counter + 1];
      int sheath_1x2v_rk3failures[sheath_1x2v_counter + 1];
      double sheath_1x2v_speciesrhs[sheath_1x2v_counter + 1];
      double sheath_1x2v_speciescollisionsrhs[sheath_1x2v_counter + 1];
      double sheath_1x2v_fieldrhs[sheath_1x2v_counter + 1];
      double sheath_1x2v_speciescollisionalmoments[sheath_1x2v_counter + 1];
      double sheath_1x2v_totalupdate[sheath_1x2v_counter + 1];
      int sheath_1x2v_memoryleakcount[sheath_1x2v_counter + 1];
      char *sheath_1x2v_memoryleaks[sheath_1x2v_counter + 1];

      for (int i = 1; i < sheath_1x2v_counter + 1; i++) {
        char *sheath_1x2v_output;
        long sheath_1x2v_file_size;
        char sheath_1x2v_buffer[128];
        snprintf(sheath_1x2v_buffer, 128, "output/rt_gk_sheath_1x2v_p1_%d.dat", i);

        FILE *sheath_1x2v_output_ptr = fopen(sheath_1x2v_buffer, "rb");
        fseek(sheath_1x2v_output_ptr, 0, SEEK_END);
        sheath_1x2v_file_size = ftell(sheath_1x2v_output_ptr);
        rewind(sheath_1x2v_output_ptr);
        sheath_1x2v_output = malloc(sheath_1x2v_file_size * (sizeof(char)));
        fread(sheath_1x2v_output, sizeof(char), sheath_1x2v_file_size, sheath_1x2v_output_ptr);
        fclose(sheath_1x2v_output_ptr);

        sheath_1x2v_updatecalls[i] = 0;
        if (strstr(sheath_1x2v_output, "Number of update calls ") != NULL) {
          char *full_substring = strstr(sheath_1x2v_output, "Number of update calls ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Number of update calls ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Number of update calls ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_1x2v_updatecalls[i] = strtol(substring, &end_ptr, 10);
        }

        sheath_1x2v_forwardeuler[i] = 0;
        if (strstr(sheath_1x2v_output, "Number of forward-Euler calls ") != NULL) {
          char *full_substring = strstr(sheath_1x2v_output, "Number of forward-Euler calls ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Number of forward-Euler calls ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Number of forward-Euler calls ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_1x2v_forwardeuler[i] = strtol(substring, &end_ptr, 10);
        }

        sheath_1x2v_rk2failures[i] = 0;
        if (strstr(sheath_1x2v_output, "Number of RK stage-2 failures ") != NULL) {
          char *full_substring = strstr(sheath_1x2v_output, "Number of RK stage-2 failures ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Number of RK stage-2 failures ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Number of RK stage-2 failures ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_1x2v_rk2failures[i] = strtol(substring, &end_ptr, 10);
        }

        sheath_1x2v_rk3failures[i] = 0;
        if (strstr(sheath_1x2v_output, "Number of RK stage-3 failures ") != NULL) {
          char *full_substring = strstr(sheath_1x2v_output, "Number of RK stage-3 failures ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Number of RK stage-3 failures ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Number of RK stage-3 failures ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_1x2v_rk3failures[i] = strtol(substring, &end_ptr, 10);
        }

        sheath_1x2v_speciesrhs[i] = 0.0;
        if (strstr(sheath_1x2v_output, "Species RHS calc took ") != NULL) {
          char *full_substring = strstr(sheath_1x2v_output, "Species RHS calc took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Species RHS calc took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Species RHS calc took ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_1x2v_speciesrhs[i] = strtod(substring, &end_ptr);
        }

        sheath_1x2v_speciescollisionsrhs[i] = 0.0;
        if (strstr(sheath_1x2v_output, "Species collisions RHS calc took ") != NULL) {
          char *full_substring = strstr(sheath_1x2v_output, "Species collisions RHS calc took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Species collisions RHS calc took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Species collisions RHS calc took ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_1x2v_speciescollisionsrhs[i] = strtod(substring, &end_ptr);
        }

        sheath_1x2v_fieldrhs[i] = 0.0;
        if (strstr(sheath_1x2v_output, "Field RHS calc took ") != NULL) {
          char *full_substring = strstr(sheath_1x2v_output, "Field RHS calc took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Field RHS calc took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Field RHS calc took ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_1x2v_fieldrhs[i] = strtod(substring, &end_ptr);
        }

        sheath_1x2v_speciescollisionalmoments[i] = 0.0;
        if (strstr(sheath_1x2v_output, "Species collisional moments took ") != NULL) {
          char *full_substring = strstr(sheath_1x2v_output, "Species collisional moments took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Species collisional moments took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Species collisional moments took ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_1x2v_speciescollisionalmoments[i] = strtod(substring, &end_ptr);
        }

        sheath_1x2v_totalupdate[i] = 0.0;
        if (strstr(sheath_1x2v_output, "Total updates took ") != NULL) {
          char *full_substring = strstr(sheath_1x2v_output, "Total updates took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Total updates took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Total updates took ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_1x2v_totalupdate[i] = strtod(substring, &end_ptr);
        }
        
        char *sheath_1x2v_temp = sheath_1x2v_output;
        sheath_1x2v_memoryleakcount[i] = 0;
        sheath_1x2v_memoryleaks[i] = (char*)malloc(1024 * sizeof(char));
        while (strstr(sheath_1x2v_temp, "0x") != NULL) {
          sheath_1x2v_temp = strstr(sheath_1x2v_temp, "0x");

          char substring[64];
          for (int j = 0; j < 64; j++) {
            substring[j] = '\0';
          }

          int substring_index = 0;
          while (sheath_1x2v_temp[substring_index] != ' ' && sheath_1x2v_temp[substring_index] != '\n') {
            substring[substring_index] = sheath_1x2v_temp[substring_index];
            substring_index += 1;
          }

          char *sheath_1x2v_temp2 = sheath_1x2v_output;
          int count = 0;
          while (strstr(sheath_1x2v_temp2, substring) != NULL) {
            sheath_1x2v_temp2 = strstr(sheath_1x2v_temp2, substring);

            count += 1;
            sheath_1x2v_temp2 += 1;
          }
          if (count == 1) {
            sheath_1x2v_memoryleakcount[i] += 1;
            sheath_1x2v_memoryleaks[i] = strcat(sheath_1x2v_memoryleaks[i], substring);
            sheath_1x2v_memoryleaks[i] = strcat(sheath_1x2v_memoryleaks[i], " ");
          }
          
          sheath_1x2v_temp += 1;
        }
      }

      for (int i = 1; i < sheath_1x2v_counter + 1; i++) {
        printf("Build number: %d\n", i);
        if (i == 1) {
          printf("Update calls: %d\n", sheath_1x2v_updatecalls[i]);
          printf("Forward-Euler calls: %d\n", sheath_1x2v_forwardeuler[i]);
          printf("RK stage-2 failures: %d\n", sheath_1x2v_rk2failures[i]);
          printf("RK stage-3 failures: %d\n", sheath_1x2v_rk3failures[i]);
          printf("Species RHS time: %f\n", sheath_1x2v_speciesrhs[i]);
          printf("Species collision RHS time: %f\n", sheath_1x2v_speciescollisionsrhs[i]);
          printf("Field RHS time: %f\n", sheath_1x2v_fieldrhs[i]);
          printf("Species collisional moments time: %f\n", sheath_1x2v_speciescollisionalmoments[i]);
          printf("Total update time: %f\n", sheath_1x2v_totalupdate[i]);
          if (sheath_1x2v_memoryleakcount[i] != 0) {
            printf("Memory leaks: " ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", sheath_1x2v_memoryleaks[i]);
          }
          else {
            printf("Memory leaks: " ANSI_COLOR_GREEN "None" ANSI_COLOR_RESET "\n");
          }
          printf("Correct: N/A\n\n");
        }
        else {
          if (sheath_1x2v_updatecalls[i] != sheath_1x2v_updatecalls[i - 1]) {
            printf("Update calls: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", sheath_1x2v_updatecalls[i]);
          }
          else {
            printf("Update calls: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", sheath_1x2v_updatecalls[i]);
          }

          if (sheath_1x2v_forwardeuler[i] != sheath_1x2v_forwardeuler[i - 1]) {
            printf("Forward-Euler calls: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", sheath_1x2v_forwardeuler[i]);
          }
          else {
            printf("Forward-Euler calls: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", sheath_1x2v_forwardeuler[i]);
          }

          if (sheath_1x2v_rk2failures[i] > sheath_1x2v_rk2failures[i - 1]) {
            printf("RK stage-2 failures: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", sheath_1x2v_rk2failures[i]);
          }
          else {
            printf("RK stage-2 failures: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", sheath_1x2v_rk2failures[i]);
          }

          if (sheath_1x2v_rk3failures[i] > sheath_1x2v_rk3failures[i - 1]) {
            printf("RK stage-3 failures: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", sheath_1x2v_rk3failures[i]);
          }
          else {
            printf("RK stage-3 failures: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", sheath_1x2v_rk3failures[i]);
          }

          if (sheath_1x2v_speciesrhs[i] > sheath_1x2v_speciesrhs[i - 1]) {
            printf("Species RHS time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", sheath_1x2v_speciesrhs[i]);
          }
          else {
            printf("Species RHS time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", sheath_1x2v_speciesrhs[i]);
          }

          if (sheath_1x2v_speciescollisionsrhs[i] > sheath_1x2v_speciescollisionsrhs[i - 1]) {
            printf("Species collision RHS time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", sheath_1x2v_speciescollisionsrhs[i]);
          }
          else {
            printf("Species collision RHS time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", sheath_1x2v_speciescollisionsrhs[i]);
          }

          if (sheath_1x2v_fieldrhs[i] > sheath_1x2v_fieldrhs[i - 1]) {
            printf("Field RHS time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", sheath_1x2v_fieldrhs[i]);
          }
          else {
            printf("Field RHS time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", sheath_1x2v_fieldrhs[i]);
          }

          if (sheath_1x2v_speciescollisionalmoments[i] > sheath_1x2v_speciescollisionalmoments[i - 1]) {
            printf("Species collisional moments time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", sheath_1x2v_speciescollisionalmoments[i]);
          }
          else {
            printf("Species collisional moments time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", sheath_1x2v_speciescollisionalmoments[i]);
          }

          if (sheath_1x2v_memoryleakcount[i] != 0) {
            printf("Memory leaks: " ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", sheath_1x2v_memoryleaks[i]);
          }
          else {
            printf("Memory leaks: " ANSI_COLOR_GREEN "None" ANSI_COLOR_RESET "\n");
          }

          if ((sheath_1x2v_updatecalls[i] != sheath_1x2v_updatecalls[i - 1]) || (sheath_1x2v_forwardeuler[i] != sheath_1x2v_forwardeuler[i - 1])
            || (sheath_1x2v_rk2failures[i] != sheath_1x2v_rk2failures[i - 1]) || (sheath_1x2v_rk3failures[i] != sheath_1x2v_rk3failures[i - 1])) {
            printf("Correct: " ANSI_COLOR_RED "No" ANSI_COLOR_RESET "\n\n");
          }
          else {
            printf("Correct: " ANSI_COLOR_GREEN "Yes" ANSI_COLOR_RESET "\n\n");
          }
        }
      }
    }
    else if (option2 == 2) {
      printf("2x2v Sheath Boundary Test with p = 1:\n\n");

      int sheath_2x2v_counter = 0;

      FILE *sheath_2x2v_counter_ptr = fopen("output/sheath_2x2v_counter.dat", "r");
      if (sheath_2x2v_counter_ptr != NULL) {
        fscanf(sheath_2x2v_counter_ptr, "%d", &sheath_2x2v_counter);
      }
      fclose(sheath_2x2v_counter_ptr);

      int sheath_2x2v_updatecalls[sheath_2x2v_counter + 1];
      int sheath_2x2v_forwardeuler[sheath_2x2v_counter + 1];
      int sheath_2x2v_rk2failures[sheath_2x2v_counter + 1];
      int sheath_2x2v_rk3failures[sheath_2x2v_counter + 1];
      double sheath_2x2v_speciesrhs[sheath_2x2v_counter + 1];
      double sheath_2x2v_speciescollisionsrhs[sheath_2x2v_counter + 1];
      double sheath_2x2v_fieldrhs[sheath_2x2v_counter + 1];
      double sheath_2x2v_speciescollisionalmoments[sheath_2x2v_counter + 1];
      double sheath_2x2v_totalupdate[sheath_2x2v_counter + 1];
      int sheath_2x2v_memoryleakcount[sheath_2x2v_counter + 1];
      char *sheath_2x2v_memoryleaks[sheath_2x2v_counter + 1];

      for (int i = 1; i < sheath_2x2v_counter + 1; i++) {
        char *sheath_2x2v_output;
        long sheath_2x2v_file_size;
        char sheath_2x2v_buffer[128];
        snprintf(sheath_2x2v_buffer, 128, "output/rt_gk_sheath_2x2v_p1_%d.dat", i);

        FILE *sheath_2x2v_output_ptr = fopen(sheath_2x2v_buffer, "rb");
        fseek(sheath_2x2v_output_ptr, 0, SEEK_END);
        sheath_2x2v_file_size = ftell(sheath_2x2v_output_ptr);
        rewind(sheath_2x2v_output_ptr);
        sheath_2x2v_output = malloc(sheath_2x2v_file_size * (sizeof(char)));
        fread(sheath_2x2v_output, sizeof(char), sheath_2x2v_file_size, sheath_2x2v_output_ptr);
        fclose(sheath_2x2v_output_ptr);

        sheath_2x2v_updatecalls[i] = 0;
        if (strstr(sheath_2x2v_output, "Number of update calls ") != NULL) {
          char *full_substring = strstr(sheath_2x2v_output, "Number of update calls ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Number of update calls ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Number of update calls ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_2x2v_updatecalls[i] = strtol(substring, &end_ptr, 10);
        }

        sheath_2x2v_forwardeuler[i] = 0;
        if (strstr(sheath_2x2v_output, "Number of forward-Euler calls ") != NULL) {
          char *full_substring = strstr(sheath_2x2v_output, "Number of forward-Euler calls ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Number of forward-Euler calls ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Number of forward-Euler calls ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_2x2v_forwardeuler[i] = strtol(substring, &end_ptr, 10);
        }

        sheath_2x2v_rk2failures[i] = 0;
        if (strstr(sheath_2x2v_output, "Number of RK stage-2 failures ") != NULL) {
          char *full_substring = strstr(sheath_2x2v_output, "Number of RK stage-2 failures ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Number of RK stage-2 failures ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Number of RK stage-2 failures ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_2x2v_rk2failures[i] = strtol(substring, &end_ptr, 10);
        }

        sheath_2x2v_rk3failures[i] = 0;
        if (strstr(sheath_2x2v_output, "Number of RK stage-3 failures ") != NULL) {
          char *full_substring = strstr(sheath_2x2v_output, "Number of RK stage-3 failures ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Number of RK stage-3 failures ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Number of RK stage-3 failures ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_2x2v_rk3failures[i] = strtol(substring, &end_ptr, 10);
        }

        sheath_2x2v_speciesrhs[i] = 0.0;
        if (strstr(sheath_2x2v_output, "Species RHS calc took ") != NULL) {
          char *full_substring = strstr(sheath_2x2v_output, "Species RHS calc took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Species RHS calc took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Species RHS calc took ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_2x2v_speciesrhs[i] = strtod(substring, &end_ptr);
        }

        sheath_2x2v_speciescollisionsrhs[i] = 0.0;
        if (strstr(sheath_2x2v_output, "Species collisions RHS calc took ") != NULL) {
          char *full_substring = strstr(sheath_2x2v_output, "Species collisions RHS calc took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Species collisions RHS calc took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Species collisions RHS calc took ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_2x2v_speciescollisionsrhs[i] = strtod(substring, &end_ptr);
        }

        sheath_2x2v_fieldrhs[i] = 0.0;
        if (strstr(sheath_2x2v_output, "Field RHS calc took ") != NULL) {
          char *full_substring = strstr(sheath_2x2v_output, "Field RHS calc took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Field RHS calc took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Field RHS calc took ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_2x2v_fieldrhs[i] = strtod(substring, &end_ptr);
        }

        sheath_2x2v_speciescollisionalmoments[i] = 0.0;
        if (strstr(sheath_2x2v_output, "Species collisional moments took ") != NULL) {
          char *full_substring = strstr(sheath_2x2v_output, "Species collisional moments took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Species collisional moments took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Species collisional moments took ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_2x2v_speciescollisionalmoments[i] = strtod(substring, &end_ptr);
        }

        sheath_2x2v_totalupdate[i] = 0.0;
        if (strstr(sheath_2x2v_output, "Total updates took ") != NULL) {
          char *full_substring = strstr(sheath_2x2v_output, "Total updates took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Total updates took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Total updates took ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_2x2v_totalupdate[i] = strtod(substring, &end_ptr);
        }

        char *sheath_2x2v_temp = sheath_2x2v_output;
        sheath_2x2v_memoryleakcount[i] = 0;
        sheath_2x2v_memoryleaks[i] = (char*)malloc(1024 * sizeof(char));
        while (strstr(sheath_2x2v_temp, "0x") != NULL) {
          sheath_2x2v_temp = strstr(sheath_2x2v_temp, "0x");

          char substring[64];
          for (int j = 0; j < 64; j++) {
            substring[j] = '\0';
          }

          int substring_index = 0;
          while (sheath_2x2v_temp[substring_index] != ' ' && sheath_2x2v_temp[substring_index] != '\n') {
            substring[substring_index] = sheath_2x2v_temp[substring_index];
            substring_index += 1;
          }

          char *sheath_2x2v_temp2 = sheath_2x2v_output;
          int count = 0;
          while (strstr(sheath_2x2v_temp2, substring) != NULL) {
            sheath_2x2v_temp2 = strstr(sheath_2x2v_temp2, substring);

            count += 1;
            sheath_2x2v_temp2 += 1;
          }
          if (count == 1) {
            sheath_2x2v_memoryleakcount[i] += 1;
            sheath_2x2v_memoryleaks[i] = strcat(sheath_2x2v_memoryleaks[i], substring);
            sheath_2x2v_memoryleaks[i] = strcat(sheath_2x2v_memoryleaks[i], " ");
          }
          
          sheath_2x2v_temp += 1;
        }
      }

      for (int i = 1; i < sheath_2x2v_counter + 1; i++) {
        printf("Build number: %d\n", i);
        if (i == 1) {
          printf("Update calls: %d\n", sheath_2x2v_updatecalls[i]);
          printf("Forward-Euler calls: %d\n", sheath_2x2v_forwardeuler[i]);
          printf("RK stage-2 failures: %d\n", sheath_2x2v_rk2failures[i]);
          printf("RK stage-3 failures: %d\n", sheath_2x2v_rk3failures[i]);
          printf("Species RHS time: %f\n", sheath_2x2v_speciesrhs[i]);
          printf("Species collision RHS time: %f\n", sheath_2x2v_speciescollisionsrhs[i]);
          printf("Field RHS time: %f\n", sheath_2x2v_fieldrhs[i]);
          printf("Species collisional moments time: %f\n", sheath_2x2v_speciescollisionalmoments[i]);
          printf("Total update time: %f\n", sheath_2x2v_totalupdate[i]);
          if (sheath_2x2v_memoryleakcount[i] != 0) {
            printf("Memory leaks: " ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", sheath_2x2v_memoryleaks[i]);
          }
          else {
            printf("Memory leaks: " ANSI_COLOR_GREEN "None" ANSI_COLOR_RESET "\n");
          }
          printf("Correct: N/A\n\n");
        }
        else {
          if (sheath_2x2v_updatecalls[i] != sheath_2x2v_updatecalls[i - 1]) {
            printf("Update calls: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", sheath_2x2v_updatecalls[i]);
          }
          else {
            printf("Update calls: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", sheath_2x2v_updatecalls[i]);
          }

          if (sheath_2x2v_forwardeuler[i] != sheath_2x2v_forwardeuler[i - 1]) {
            printf("Forward-Euler calls: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", sheath_2x2v_forwardeuler[i]);
          }
          else {
            printf("Forward-Euler calls: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", sheath_2x2v_forwardeuler[i]);
          }

          if (sheath_2x2v_rk2failures[i] > sheath_2x2v_rk2failures[i - 1]) {
            printf("RK stage-2 failures: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", sheath_2x2v_rk2failures[i]);
          }
          else {
            printf("RK stage-2 failures: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", sheath_2x2v_rk2failures[i]);
          }

          if (sheath_2x2v_rk3failures[i] > sheath_2x2v_rk3failures[i - 1]) {
            printf("RK stage-3 failures: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", sheath_2x2v_rk3failures[i]);
          }
          else {
            printf("RK stage-3 failures: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", sheath_2x2v_rk3failures[i]);
          }

          if (sheath_2x2v_speciesrhs[i] > sheath_2x2v_speciesrhs[i - 1]) {
            printf("Species RHS time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", sheath_2x2v_speciesrhs[i]);
          }
          else {
            printf("Species RHS time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", sheath_2x2v_speciesrhs[i]);
          }

          if (sheath_2x2v_speciescollisionsrhs[i] > sheath_2x2v_speciescollisionsrhs[i - 1]) {
            printf("Species collision RHS time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", sheath_2x2v_speciescollisionsrhs[i]);
          }
          else {
            printf("Species collision RHS time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", sheath_2x2v_speciescollisionsrhs[i]);
          }

          if (sheath_2x2v_fieldrhs[i] > sheath_2x2v_fieldrhs[i - 1]) {
            printf("Field RHS time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", sheath_2x2v_fieldrhs[i]);
          }
          else {
            printf("Field RHS time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", sheath_2x2v_fieldrhs[i]);
          }

          if (sheath_2x2v_speciescollisionalmoments[i] > sheath_2x2v_speciescollisionalmoments[i - 1]) {
            printf("Species collisional moments time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", sheath_2x2v_speciescollisionalmoments[i]);
          }
          else {
            printf("Species collisional moments time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", sheath_2x2v_speciescollisionalmoments[i]);
          }

          if (sheath_2x2v_totalupdate[i] > sheath_2x2v_totalupdate[i - 1]) {
            printf("Total update time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", sheath_2x2v_totalupdate[i]);
          }
          else {
            printf("Total update time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", sheath_2x2v_totalupdate[i]);
          }

          if (sheath_2x2v_memoryleakcount[i] != 0) {
            printf("Memory leaks: " ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", sheath_2x2v_memoryleaks[i]);
          }
          else {
            printf("Memory leaks: " ANSI_COLOR_GREEN "None" ANSI_COLOR_RESET "\n");
          }

          if ((sheath_2x2v_updatecalls[i] != sheath_2x2v_updatecalls[i - 1]) || (sheath_2x2v_forwardeuler[i] != sheath_2x2v_forwardeuler[i - 1])
            || (sheath_2x2v_rk2failures[i] != sheath_2x2v_rk2failures[i - 1]) || (sheath_2x2v_rk3failures[i] != sheath_2x2v_rk3failures[i - 1])) {
            printf("Correct: " ANSI_COLOR_RED "No" ANSI_COLOR_RESET "\n\n");
          }
          else {
            printf("Correct: " ANSI_COLOR_GREEN "Yes" ANSI_COLOR_RESET "\n\n");
          }
        }
      }
    }
    else if (option2 == 3) {
      printf("3x2v Sheath Boundary Test with p = 1:\n\n");

      int sheath_3x2v_counter = 0;

      FILE *sheath_3x2v_counter_ptr = fopen("output/sheath_3x2v_counter.dat", "r");
      if (sheath_3x2v_counter_ptr != NULL) {
        fscanf(sheath_3x2v_counter_ptr, "%d", &sheath_3x2v_counter);
      }
      fclose(sheath_3x2v_counter_ptr);

      int sheath_3x2v_updatecalls[sheath_3x2v_counter + 1];
      int sheath_3x2v_forwardeuler[sheath_3x2v_counter + 1];
      int sheath_3x2v_rk2failures[sheath_3x2v_counter + 1];
      int sheath_3x2v_rk3failures[sheath_3x2v_counter + 1];
      double sheath_3x2v_speciesrhs[sheath_3x2v_counter + 1];
      double sheath_3x2v_speciescollisionsrhs[sheath_3x2v_counter + 1];
      double sheath_3x2v_fieldrhs[sheath_3x2v_counter + 1];
      double sheath_3x2v_speciescollisionalmoments[sheath_3x2v_counter + 1];
      double sheath_3x2v_totalupdate[sheath_3x2v_counter + 1];
      int sheath_3x2v_memoryleakcount[sheath_3x2v_counter + 1];
      char *sheath_3x2v_memoryleaks[sheath_3x2v_counter + 1];

      for (int i = 1; i < sheath_3x2v_counter + 1; i++) {
        char *sheath_3x2v_output;
        long sheath_3x2v_file_size;
        char sheath_3x2v_buffer[128];
        snprintf(sheath_3x2v_buffer, 128, "output/rt_gk_sheath_3x2v_p1_%d.dat", i);

        FILE *sheath_3x2v_output_ptr = fopen(sheath_3x2v_buffer, "rb");
        fseek(sheath_3x2v_output_ptr, 0, SEEK_END);
        sheath_3x2v_file_size = ftell(sheath_3x2v_output_ptr);
        rewind(sheath_3x2v_output_ptr);
        sheath_3x2v_output = malloc(sheath_3x2v_file_size * (sizeof(char)));
        fread(sheath_3x2v_output, sizeof(char), sheath_3x2v_file_size, sheath_3x2v_output_ptr);
        fclose(sheath_3x2v_output_ptr);

        sheath_3x2v_updatecalls[i] = 0;
        if (strstr(sheath_3x2v_output, "Number of update calls ") != NULL) {
          char *full_substring = strstr(sheath_3x2v_output, "Number of update calls ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Number of update calls ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Number of update calls ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_3x2v_updatecalls[i] = strtol(substring, &end_ptr, 10);
        }

        sheath_3x2v_forwardeuler[i] = 0;
        if (strstr(sheath_3x2v_output, "Number of forward-Euler calls ") != NULL) {
          char *full_substring = strstr(sheath_3x2v_output, "Number of forward-Euler calls ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Number of forward-Euler calls ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Number of forward-Euler calls ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_3x2v_forwardeuler[i] = strtol(substring, &end_ptr, 10);
        }

        sheath_3x2v_rk2failures[i] = 0;
        if (strstr(sheath_3x2v_output, "Number of RK stage-2 failures ") != NULL) {
          char *full_substring = strstr(sheath_3x2v_output, "Number of RK stage-2 failures ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Number of RK stage-2 failures ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Number of RK stage-2 failures ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_3x2v_rk2failures[i] = strtol(substring, &end_ptr, 10);
        }

        sheath_3x2v_rk3failures[i] = 0;
        if (strstr(sheath_3x2v_output, "Number of RK stage-3 failures ") != NULL) {
          char *full_substring = strstr(sheath_3x2v_output, "Number of RK stage-3 failures ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Number of RK stage-3 failures ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Number of RK stage-3 failures ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_3x2v_rk3failures[i] = strtol(substring, &end_ptr, 10);
        }

        sheath_3x2v_speciesrhs[i] = 0.0;
        if (strstr(sheath_3x2v_output, "Species RHS calc took ") != NULL) {
          char *full_substring = strstr(sheath_3x2v_output, "Species RHS calc took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Species RHS calc took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Species RHS calc took ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_3x2v_speciesrhs[i] = strtod(substring, &end_ptr);
        }

        sheath_3x2v_speciescollisionsrhs[i] = 0.0;
        if (strstr(sheath_3x2v_output, "Species collisions RHS calc took ") != NULL) {
          char *full_substring = strstr(sheath_3x2v_output, "Species collisions RHS calc took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Species collisions RHS calc took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Species collisions RHS calc took ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_3x2v_speciescollisionsrhs[i] = strtod(substring, &end_ptr);
        }

        sheath_3x2v_fieldrhs[i] = 0.0;
        if (strstr(sheath_3x2v_output, "Field RHS calc took ") != NULL) {
          char *full_substring = strstr(sheath_3x2v_output, "Field RHS calc took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Field RHS calc took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Field RHS calc took ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_3x2v_fieldrhs[i] = strtod(substring, &end_ptr);
        }

        sheath_3x2v_speciescollisionalmoments[i] = 0.0;
        if (strstr(sheath_3x2v_output, "Species collisional moments took ") != NULL) {
          char *full_substring = strstr(sheath_3x2v_output, "Species collisional moments took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Species collisional moments took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Species collisional moments took ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_3x2v_speciescollisionalmoments[i] = strtod(substring, &end_ptr);
        }

        sheath_3x2v_totalupdate[i] = 0.0;
        if (strstr(sheath_3x2v_output, "Total updates took ") != NULL) {
          char *full_substring = strstr(sheath_3x2v_output, "Total updates took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Total updates took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Total updates took ")];
            substring_index += 1;
          }

          char *end_ptr;
          sheath_3x2v_totalupdate[i] = strtod(substring, &end_ptr);
        }

        char *sheath_3x2v_temp = sheath_3x2v_output;
        sheath_3x2v_memoryleakcount[i] = 0;
        sheath_3x2v_memoryleaks[i] = (char*)malloc(1024 * sizeof(char));
        while (strstr(sheath_3x2v_temp, "0x") != NULL) {
          sheath_3x2v_temp = strstr(sheath_3x2v_temp, "0x");

          char substring[64];
          for (int j = 0; j < 64; j++) {
            substring[j] = '\0';
          }

          int substring_index = 0;
          while (sheath_3x2v_temp[substring_index] != ' ' && sheath_3x2v_temp[substring_index] != '\n') {
            substring[substring_index] = sheath_3x2v_temp[substring_index];
            substring_index += 1;
          }

          char *sheath_3x2v_temp2 = sheath_3x2v_output;
          int count = 0;
          while (strstr(sheath_3x2v_temp2, substring) != NULL) {
            sheath_3x2v_temp2 = strstr(sheath_3x2v_temp2, substring);

            count += 1;
            sheath_3x2v_temp2 += 1;
          }
          if (count == 1) {
            sheath_3x2v_memoryleakcount[i] += 1;
            sheath_3x2v_memoryleaks[i] = strcat(sheath_3x2v_memoryleaks[i], substring);
            sheath_3x2v_memoryleaks[i] = strcat(sheath_3x2v_memoryleaks[i], " ");
          }
          
          sheath_3x2v_temp += 1;
        }
      }

      for (int i = 1; i < sheath_3x2v_counter + 1; i++) {
        printf("Build number: %d\n", i);
        if (i == 1) {
          printf("Update calls: %d\n", sheath_3x2v_updatecalls[i]);
          printf("Forward-Euler calls: %d\n", sheath_3x2v_forwardeuler[i]);
          printf("RK stage-2 failures: %d\n", sheath_3x2v_rk2failures[i]);
          printf("RK stage-3 failures: %d\n", sheath_3x2v_rk3failures[i]);
          printf("Species RHS time: %f\n", sheath_3x2v_speciesrhs[i]);
          printf("Species collision RHS time: %f\n", sheath_3x2v_speciescollisionsrhs[i]);
          printf("Field RHS time: %f\n", sheath_3x2v_fieldrhs[i]);
          printf("Species collisional moments time: %f\n", sheath_3x2v_speciescollisionalmoments[i]);
          printf("Total update time: %f\n", sheath_3x2v_totalupdate[i]);
          if (sheath_3x2v_memoryleakcount[i] != 0) {
            printf("Memory leaks: " ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", sheath_3x2v_memoryleaks[i]);
          }
          else {
            printf("Memory leaks: " ANSI_COLOR_GREEN "None" ANSI_COLOR_RESET "\n");
          }
          printf("Correct: N/A\n\n");
        }
        else {
          if (sheath_3x2v_updatecalls[i] != sheath_3x2v_updatecalls[i - 1]) {
            printf("Update calls: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", sheath_3x2v_updatecalls[i]);
          }
          else {
            printf("Update calls: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", sheath_3x2v_updatecalls[i]);
          }

          if (sheath_3x2v_forwardeuler[i] != sheath_3x2v_forwardeuler[i - 1]) {
            printf("Forward-Euler calls: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", sheath_3x2v_forwardeuler[i]);
          }
          else {
            printf("Forward-Euler calls: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", sheath_3x2v_forwardeuler[i]);
          }

          if (sheath_3x2v_rk2failures[i] > sheath_3x2v_rk2failures[i - 1]) {
            printf("RK stage-2 failures: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", sheath_3x2v_rk2failures[i]);
          }
          else {
            printf("RK stage-2 failures: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", sheath_3x2v_rk2failures[i]);
          }

          if (sheath_3x2v_rk3failures[i] > sheath_3x2v_rk3failures[i - 1]) {
            printf("RK stage-3 failures: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", sheath_3x2v_rk3failures[i]);
          }
          else {
            printf("RK stage-3 failures: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", sheath_3x2v_rk3failures[i]);
          }

          if (sheath_3x2v_speciesrhs[i] > sheath_3x2v_speciesrhs[i - 1]) {
            printf("Species RHS time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", sheath_3x2v_speciesrhs[i]);
          }
          else {
            printf("Species RHS time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", sheath_3x2v_speciesrhs[i]);
          }

          if (sheath_3x2v_speciescollisionsrhs[i] > sheath_3x2v_speciescollisionsrhs[i - 1]) {
            printf("Species collision RHS time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", sheath_3x2v_speciescollisionsrhs[i]);
          }
          else {
            printf("Species collision RHS time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", sheath_3x2v_speciescollisionsrhs[i]);
          }

          if (sheath_3x2v_fieldrhs[i] > sheath_3x2v_fieldrhs[i - 1]) {
            printf("Field RHS time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", sheath_3x2v_fieldrhs[i]);
          }
          else {
            printf("Field RHS time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", sheath_3x2v_fieldrhs[i]);
          }

          if (sheath_3x2v_speciescollisionalmoments[i] > sheath_3x2v_speciescollisionalmoments[i - 1]) {
            printf("Species collisional moments time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", sheath_3x2v_speciescollisionalmoments[i]);
          }
          else {
            printf("Species collisional moments time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", sheath_3x2v_speciescollisionalmoments[i]);
          }

          if (sheath_3x2v_totalupdate[i] > sheath_3x2v_totalupdate[i - 1]) {
            printf("Total update time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", sheath_3x2v_totalupdate[i]);
          }
          else {
            printf("Total update time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", sheath_3x2v_totalupdate[i]);
          }

          if (sheath_3x2v_memoryleakcount[i] != 0) {
            printf("Memory leaks: " ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", sheath_3x2v_memoryleaks[i]);
          }
          else {
            printf("Memory leaks: " ANSI_COLOR_GREEN "None" ANSI_COLOR_RESET "\n");
          }

          if ((sheath_3x2v_updatecalls[i] != sheath_3x2v_updatecalls[i - 1]) || (sheath_3x2v_forwardeuler[i] != sheath_3x2v_forwardeuler[i - 1])
            || (sheath_3x2v_rk2failures[i] != sheath_3x2v_rk2failures[i - 1]) || (sheath_3x2v_rk3failures[i] != sheath_3x2v_rk3failures[i - 1])) {
            printf("Correct: " ANSI_COLOR_RED "No" ANSI_COLOR_RESET "\n\n");
          }
          else {
            printf("Correct: " ANSI_COLOR_GREEN "Yes" ANSI_COLOR_RESET "\n\n");
          }
        }
      }
    }
    else if (option2 == 4) {
      printf("3x2v LAPD Test (in Cartesian coordinates) with p = 1:\n\n");

      int lapd_cart_3x2v_counter = 0;

      FILE *lapd_cart_3x2v_counter_ptr = fopen("output/lapd_cart_3x2v_counter.dat", "r");
      if (lapd_cart_3x2v_counter_ptr != NULL) {
        fscanf(lapd_cart_3x2v_counter_ptr, "%d", &lapd_cart_3x2v_counter);
      }
      fclose(lapd_cart_3x2v_counter_ptr);

      int lapd_cart_3x2v_updatecalls[lapd_cart_3x2v_counter + 1];
      int lapd_cart_3x2v_forwardeuler[lapd_cart_3x2v_counter + 1];
      int lapd_cart_3x2v_rk2failures[lapd_cart_3x2v_counter + 1];
      int lapd_cart_3x2v_rk3failures[lapd_cart_3x2v_counter + 1];
      double lapd_cart_3x2v_speciesrhs[lapd_cart_3x2v_counter + 1];
      double lapd_cart_3x2v_speciescollisionsrhs[lapd_cart_3x2v_counter + 1];
      double lapd_cart_3x2v_fieldrhs[lapd_cart_3x2v_counter + 1];
      double lapd_cart_3x2v_speciescollisionalmoments[lapd_cart_3x2v_counter + 1];
      double lapd_cart_3x2v_totalupdate[lapd_cart_3x2v_counter + 1];
      int lapd_cart_3x2v_memoryleakcount[lapd_cart_3x2v_counter + 1];
      char *lapd_cart_3x2v_memoryleaks[lapd_cart_3x2v_counter + 1];

      for (int i = 1; i < lapd_cart_3x2v_counter + 1; i++) {
        char *lapd_cart_3x2v_output;
        long lapd_cart_3x2v_file_size;
        char lapd_cart_3x2v_buffer[128];
        snprintf(lapd_cart_3x2v_buffer, 128, "output/rt_gk_lapd_cart_3x2v_p1_%d.dat", i);

        FILE *lapd_cart_3x2v_output_ptr = fopen(lapd_cart_3x2v_buffer, "rb");
        fseek(lapd_cart_3x2v_output_ptr, 0, SEEK_END);
        lapd_cart_3x2v_file_size = ftell(lapd_cart_3x2v_output_ptr);
        rewind(lapd_cart_3x2v_output_ptr);
        lapd_cart_3x2v_output = malloc(lapd_cart_3x2v_file_size * (sizeof(char)));
        fread(lapd_cart_3x2v_output, sizeof(char), lapd_cart_3x2v_file_size, lapd_cart_3x2v_output_ptr);
        fclose(lapd_cart_3x2v_output_ptr);

        lapd_cart_3x2v_updatecalls[i] = 0;
        if (strstr(lapd_cart_3x2v_output, "Number of update calls ") != NULL) {
          char *full_substring = strstr(lapd_cart_3x2v_output, "Number of update calls ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Number of update calls ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Number of update calls ")];
            substring_index += 1;
          }

          char *end_ptr;
          lapd_cart_3x2v_updatecalls[i] = strtol(substring, &end_ptr, 10);
        }

        lapd_cart_3x2v_forwardeuler[i] = 0;
        if (strstr(lapd_cart_3x2v_output, "Number of forward-Euler calls ") != NULL) {
          char *full_substring = strstr(lapd_cart_3x2v_output, "Number of forward-Euler calls ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Number of forward-Euler calls ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Number of forward-Euler calls ")];
            substring_index += 1;
          }

          char *end_ptr;
          lapd_cart_3x2v_forwardeuler[i] = strtol(substring, &end_ptr, 10);
        }

        lapd_cart_3x2v_rk2failures[i] = 0;
        if (strstr(lapd_cart_3x2v_output, "Number of RK stage-2 failures ") != NULL) {
          char *full_substring = strstr(lapd_cart_3x2v_output, "Number of RK stage-2 failures ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Number of RK stage-2 failures ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Number of RK stage-2 failures ")];
            substring_index += 1;
          }

          char *end_ptr;
          lapd_cart_3x2v_rk2failures[i] = strtol(substring, &end_ptr, 10);
        }

        lapd_cart_3x2v_rk3failures[i] = 0;
        if (strstr(lapd_cart_3x2v_output, "Number of RK stage-3 failures ") != NULL) {
          char *full_substring = strstr(lapd_cart_3x2v_output, "Number of RK stage-3 failures ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Number of RK stage-3 failures ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Number of RK stage-3 failures ")];
            substring_index += 1;
          }

          char *end_ptr;
          lapd_cart_3x2v_rk3failures[i] = strtol(substring, &end_ptr, 10);
        }

        lapd_cart_3x2v_speciesrhs[i] = 0.0;
        if (strstr(lapd_cart_3x2v_output, "Species RHS calc took ") != NULL) {
          char *full_substring = strstr(lapd_cart_3x2v_output, "Species RHS calc took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Species RHS calc took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Species RHS calc took ")];
            substring_index += 1;
          }

          char *end_ptr;
          lapd_cart_3x2v_speciesrhs[i] = strtod(substring, &end_ptr);
        }

        lapd_cart_3x2v_speciescollisionsrhs[i] = 0.0;
        if (strstr(lapd_cart_3x2v_output, "Species collisions RHS calc took ") != NULL) {
          char *full_substring = strstr(lapd_cart_3x2v_output, "Species collisions RHS calc took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Species collisions RHS calc took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Species collisions RHS calc took ")];
            substring_index += 1;
          }

          char *end_ptr;
          lapd_cart_3x2v_speciescollisionsrhs[i] = strtod(substring, &end_ptr);
        }

        lapd_cart_3x2v_fieldrhs[i] = 0.0;
        if (strstr(lapd_cart_3x2v_output, "Field RHS calc took ") != NULL) {
          char *full_substring = strstr(lapd_cart_3x2v_output, "Field RHS calc took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Field RHS calc took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Field RHS calc took ")];
            substring_index += 1;
          }

          char *end_ptr;
          lapd_cart_3x2v_fieldrhs[i] = strtod(substring, &end_ptr);
        }

        lapd_cart_3x2v_speciescollisionalmoments[i] = 0.0;
        if (strstr(lapd_cart_3x2v_output, "Species collisional moments took ") != NULL) {
          char *full_substring = strstr(lapd_cart_3x2v_output, "Species collisional moments took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Species collisional moments took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Species collisional moments took ")];
            substring_index += 1;
          }

          char *end_ptr;
          lapd_cart_3x2v_speciescollisionalmoments[i] = strtod(substring, &end_ptr);
        }

        lapd_cart_3x2v_totalupdate[i] = 0.0;
        if (strstr(lapd_cart_3x2v_output, "Total updates took ") != NULL) {
          char *full_substring = strstr(lapd_cart_3x2v_output, "Total updates took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Total updates took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Total updates took ")];
            substring_index += 1;
          }

          char *end_ptr;
          lapd_cart_3x2v_totalupdate[i] = strtod(substring, &end_ptr);
        }

        char *lapd_cart_3x2v_temp = lapd_cart_3x2v_output;
        lapd_cart_3x2v_memoryleakcount[i] = 0;
        lapd_cart_3x2v_memoryleaks[i] = (char*)malloc(1024 * sizeof(char));
        while (strstr(lapd_cart_3x2v_temp, "0x") != NULL) {
          lapd_cart_3x2v_temp = strstr(lapd_cart_3x2v_temp, "0x");

          char substring[64];
          for (int j = 0; j < 64; j++) {
            substring[j] = '\0';
          }

          int substring_index = 0;
          while (lapd_cart_3x2v_temp[substring_index] != ' ' && lapd_cart_3x2v_temp[substring_index] != '\n') {
            substring[substring_index] = lapd_cart_3x2v_temp[substring_index];
            substring_index += 1;
          }

          char *lapd_cart_3x2v_temp2 = lapd_cart_3x2v_output;
          int count = 0;
          while (strstr(lapd_cart_3x2v_temp2, substring) != NULL) {
            lapd_cart_3x2v_temp2 = strstr(lapd_cart_3x2v_temp2, substring);

            count += 1;
            lapd_cart_3x2v_temp2 += 1;
          }
          if (count == 1) {
            lapd_cart_3x2v_memoryleakcount[i] += 1;
            lapd_cart_3x2v_memoryleaks[i] = strcat(lapd_cart_3x2v_memoryleaks[i], substring);
            lapd_cart_3x2v_memoryleaks[i] = strcat(lapd_cart_3x2v_memoryleaks[i], " ");
          }
          
          lapd_cart_3x2v_temp += 1;
        }
      }

      for (int i = 1; i < lapd_cart_3x2v_counter + 1; i++) {
        printf("Build number: %d\n", i);
        if (i == 1) {
          printf("Update calls: %d\n", lapd_cart_3x2v_updatecalls[i]);
          printf("Forward-Euler calls: %d\n", lapd_cart_3x2v_forwardeuler[i]);
          printf("RK stage-2 failures: %d\n", lapd_cart_3x2v_rk2failures[i]);
          printf("RK stage-3 failures: %d\n", lapd_cart_3x2v_rk3failures[i]);
          printf("Species RHS time: %f\n", lapd_cart_3x2v_speciesrhs[i]);
          printf("Species collision RHS time: %f\n", lapd_cart_3x2v_speciescollisionsrhs[i]);
          printf("Field RHS time: %f\n", lapd_cart_3x2v_fieldrhs[i]);
          printf("Species collisional moments time: %f\n", lapd_cart_3x2v_speciescollisionalmoments[i]);
          printf("Total update time: %f\n", lapd_cart_3x2v_totalupdate[i]);
          if (lapd_cart_3x2v_memoryleakcount[i] != 0) {
            printf("Memory leaks: " ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", lapd_cart_3x2v_memoryleaks[i]);
          }
          else {
            printf("Memory leaks: " ANSI_COLOR_GREEN "None" ANSI_COLOR_RESET "\n");
          }
          printf("Correct: N/A\n\n");
        }
        else {
          if (lapd_cart_3x2v_updatecalls[i] != lapd_cart_3x2v_updatecalls[i - 1]) {
            printf("Update calls: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", lapd_cart_3x2v_updatecalls[i]);
          }
          else {
            printf("Update calls: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", lapd_cart_3x2v_updatecalls[i]);
          }

          if (lapd_cart_3x2v_forwardeuler[i] != lapd_cart_3x2v_forwardeuler[i - 1]) {
            printf("Forward-Euler calls: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", lapd_cart_3x2v_forwardeuler[i]);
          }
          else {
            printf("Forward-Euler calls: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", lapd_cart_3x2v_forwardeuler[i]);
          }

          if (lapd_cart_3x2v_rk2failures[i] > lapd_cart_3x2v_rk2failures[i - 1]) {
            printf("RK stage-2 failures: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", lapd_cart_3x2v_rk2failures[i]);
          }
          else {
            printf("RK stage-2 failures: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", lapd_cart_3x2v_rk2failures[i]);
          }

          if (lapd_cart_3x2v_rk3failures[i] > lapd_cart_3x2v_rk3failures[i - 1]) {
            printf("RK stage-3 failures: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", lapd_cart_3x2v_rk3failures[i]);
          }
          else {
            printf("RK stage-3 failures: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", lapd_cart_3x2v_rk3failures[i]);
          }

          if (lapd_cart_3x2v_speciesrhs[i] > lapd_cart_3x2v_speciesrhs[i - 1]) {
            printf("Species RHS time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", lapd_cart_3x2v_speciesrhs[i]);
          }
          else {
            printf("Species RHS time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", lapd_cart_3x2v_speciesrhs[i]);
          }

          if (lapd_cart_3x2v_speciescollisionsrhs[i] > lapd_cart_3x2v_speciescollisionsrhs[i - 1]) {
            printf("Species collision RHS time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", lapd_cart_3x2v_speciescollisionsrhs[i]);
          }
          else {
            printf("Species collision RHS time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", lapd_cart_3x2v_speciescollisionsrhs[i]);
          }

          if (lapd_cart_3x2v_fieldrhs[i] > lapd_cart_3x2v_fieldrhs[i - 1]) {
            printf("Field RHS time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", lapd_cart_3x2v_fieldrhs[i]);
          }
          else {
            printf("Field RHS time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", lapd_cart_3x2v_fieldrhs[i]);
          }

          if (lapd_cart_3x2v_speciescollisionalmoments[i] > lapd_cart_3x2v_speciescollisionalmoments[i - 1]) {
            printf("Species collisional moments time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", lapd_cart_3x2v_speciescollisionalmoments[i]);
          }
          else {
            printf("Species collisional moments time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", lapd_cart_3x2v_speciescollisionalmoments[i]);
          }

          if (lapd_cart_3x2v_totalupdate[i] > lapd_cart_3x2v_totalupdate[i - 1]) {
            printf("Total update time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", lapd_cart_3x2v_totalupdate[i]);
          }
          else {
            printf("Total update time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", lapd_cart_3x2v_totalupdate[i]);
          }

          if (lapd_cart_3x2v_memoryleakcount[i] != 0) {
            printf("Memory leaks: " ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", lapd_cart_3x2v_memoryleaks[i]);
          }
          else {
            printf("Memory leaks: " ANSI_COLOR_GREEN "None" ANSI_COLOR_RESET "\n");
          }

          if ((lapd_cart_3x2v_updatecalls[i] != lapd_cart_3x2v_updatecalls[i - 1]) || (lapd_cart_3x2v_forwardeuler[i] != lapd_cart_3x2v_forwardeuler[i - 1])
            || (lapd_cart_3x2v_rk2failures[i] != lapd_cart_3x2v_rk2failures[i - 1]) || (lapd_cart_3x2v_rk3failures[i] != lapd_cart_3x2v_rk3failures[i - 1])) {
            printf("Correct: " ANSI_COLOR_RED "No" ANSI_COLOR_RESET "\n\n");
          }
          else {
            printf("Correct: " ANSI_COLOR_GREEN "Yes" ANSI_COLOR_RESET "\n\n");
          }
        }
      }
    }
    else if (option2 == 5) {
      printf("3x2v LAPD Test (in cylindrical coordinates) with p = 1:\n\n");

      int lapd_cyl_3x2v_counter = 0;

      FILE *lapd_cyl_3x2v_counter_ptr = fopen("output/lapd_cyl_3x2v_counter.dat", "r");
      if (lapd_cyl_3x2v_counter_ptr != NULL) {
        fscanf(lapd_cyl_3x2v_counter_ptr, "%d", &lapd_cyl_3x2v_counter);
      }
      fclose(lapd_cyl_3x2v_counter_ptr);

      int lapd_cyl_3x2v_updatecalls[lapd_cyl_3x2v_counter + 1];
      int lapd_cyl_3x2v_forwardeuler[lapd_cyl_3x2v_counter + 1];
      int lapd_cyl_3x2v_rk2failures[lapd_cyl_3x2v_counter + 1];
      int lapd_cyl_3x2v_rk3failures[lapd_cyl_3x2v_counter + 1];
      double lapd_cyl_3x2v_speciesrhs[lapd_cyl_3x2v_counter + 1];
      double lapd_cyl_3x2v_speciescollisionsrhs[lapd_cyl_3x2v_counter + 1];
      double lapd_cyl_3x2v_fieldrhs[lapd_cyl_3x2v_counter + 1];
      double lapd_cyl_3x2v_speciescollisionalmoments[lapd_cyl_3x2v_counter + 1];
      double lapd_cyl_3x2v_totalupdate[lapd_cyl_3x2v_counter + 1];
      int lapd_cyl_3x2v_memoryleakcount[lapd_cyl_3x2v_counter + 1];
      char *lapd_cyl_3x2v_memoryleaks[lapd_cyl_3x2v_counter + 1];

      for (int i = 1; i < lapd_cyl_3x2v_counter + 1; i++) {
        char *lapd_cyl_3x2v_output;
        long lapd_cyl_3x2v_file_size;
        char lapd_cyl_3x2v_buffer[128];
        snprintf(lapd_cyl_3x2v_buffer, 128, "output/rt_gk_lapd_cyl_3x2v_p1_%d.dat", i);

        FILE *lapd_cyl_3x2v_output_ptr = fopen(lapd_cyl_3x2v_buffer, "rb");
        fseek(lapd_cyl_3x2v_output_ptr, 0, SEEK_END);
        lapd_cyl_3x2v_file_size = ftell(lapd_cyl_3x2v_output_ptr);
        rewind(lapd_cyl_3x2v_output_ptr);
        lapd_cyl_3x2v_output = malloc(lapd_cyl_3x2v_file_size * (sizeof(char)));
        fread(lapd_cyl_3x2v_output, sizeof(char), lapd_cyl_3x2v_file_size, lapd_cyl_3x2v_output_ptr);
        fclose(lapd_cyl_3x2v_output_ptr);

        lapd_cyl_3x2v_updatecalls[i] = 0;
        if (strstr(lapd_cyl_3x2v_output, "Number of update calls ") != NULL) {
          char *full_substring = strstr(lapd_cyl_3x2v_output, "Number of update calls ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Number of update calls ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Number of update calls ")];
            substring_index += 1;
          }

          char *end_ptr;
          lapd_cyl_3x2v_updatecalls[i] = strtol(substring, &end_ptr, 10);
        }

        lapd_cyl_3x2v_forwardeuler[i] = 0;
        if (strstr(lapd_cyl_3x2v_output, "Number of forward-Euler calls ") != NULL) {
          char *full_substring = strstr(lapd_cyl_3x2v_output, "Number of forward-Euler calls ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Number of forward-Euler calls ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Number of forward-Euler calls ")];
            substring_index += 1;
          }

          char *end_ptr;
          lapd_cyl_3x2v_forwardeuler[i] = strtol(substring, &end_ptr, 10);
        }

        lapd_cyl_3x2v_rk2failures[i] = 0;
        if (strstr(lapd_cyl_3x2v_output, "Number of RK stage-2 failures ") != NULL) {
          char *full_substring = strstr(lapd_cyl_3x2v_output, "Number of RK stage-2 failures ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Number of RK stage-2 failures ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Number of RK stage-2 failures ")];
            substring_index += 1;
          }

          char *end_ptr;
          lapd_cyl_3x2v_rk2failures[i] = strtol(substring, &end_ptr, 10);
        }

        lapd_cyl_3x2v_rk3failures[i] = 0;
        if (strstr(lapd_cyl_3x2v_output, "Number of RK stage-3 failures ") != NULL) {
          char *full_substring = strstr(lapd_cyl_3x2v_output, "Number of RK stage-3 failures ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Number of RK stage-3 failures ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Number of RK stage-3 failures ")];
            substring_index += 1;
          }

          char *end_ptr;
          lapd_cyl_3x2v_rk3failures[i] = strtol(substring, &end_ptr, 10);
        }

        lapd_cyl_3x2v_speciesrhs[i] = 0.0;
        if (strstr(lapd_cyl_3x2v_output, "Species RHS calc took ") != NULL) {
          char *full_substring = strstr(lapd_cyl_3x2v_output, "Species RHS calc took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Species RHS calc took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Species RHS calc took ")];
            substring_index += 1;
          }

          char *end_ptr;
          lapd_cyl_3x2v_speciesrhs[i] = strtod(substring, &end_ptr);
        }

        lapd_cyl_3x2v_speciescollisionsrhs[i] = 0.0;
        if (strstr(lapd_cyl_3x2v_output, "Species collisions RHS calc took ") != NULL) {
          char *full_substring = strstr(lapd_cyl_3x2v_output, "Species collisions RHS calc took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Species collisions RHS calc took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Species collisions RHS calc took ")];
            substring_index += 1;
          }

          char *end_ptr;
          lapd_cyl_3x2v_speciescollisionsrhs[i] = strtod(substring, &end_ptr);
        }

        lapd_cyl_3x2v_fieldrhs[i] = 0.0;
        if (strstr(lapd_cyl_3x2v_output, "Field RHS calc took ") != NULL) {
          char *full_substring = strstr(lapd_cyl_3x2v_output, "Field RHS calc took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Field RHS calc took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Field RHS calc took ")];
            substring_index += 1;
          }

          char *end_ptr;
          lapd_cyl_3x2v_fieldrhs[i] = strtod(substring, &end_ptr);
        }

        lapd_cyl_3x2v_speciescollisionalmoments[i] = 0.0;
        if (strstr(lapd_cyl_3x2v_output, "Species collisional moments took ") != NULL) {
          char *full_substring = strstr(lapd_cyl_3x2v_output, "Species collisional moments took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Species collisional moments took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Species collisional moments took ")];
            substring_index += 1;
          }

          char *end_ptr;
          lapd_cyl_3x2v_speciescollisionalmoments[i] = strtod(substring, &end_ptr);
        }

        lapd_cyl_3x2v_totalupdate[i] = 0.0;
        if (strstr(lapd_cyl_3x2v_output, "Total updates took ") != NULL) {
          char *full_substring = strstr(lapd_cyl_3x2v_output, "Total updates took ");
          char substring[64];
          int substring_index = 0;

          while (full_substring[substring_index + strlen("Total updates took ")] != '\n') {
            substring[substring_index] = full_substring[substring_index + strlen("Total updates took ")];
            substring_index += 1;
          }

          char *end_ptr;
          lapd_cyl_3x2v_totalupdate[i] = strtod(substring, &end_ptr);
        }

        char *lapd_cyl_3x2v_temp = lapd_cyl_3x2v_output;
        lapd_cyl_3x2v_memoryleakcount[i] = 0;
        lapd_cyl_3x2v_memoryleaks[i] = (char*)malloc(1024 * sizeof(char));
        while (strstr(lapd_cyl_3x2v_temp, "0x") != NULL) {
          lapd_cyl_3x2v_temp = strstr(lapd_cyl_3x2v_temp, "0x");

          char substring[64];
          for (int j = 0; j < 64; j++) {
            substring[j] = '\0';
          }

          int substring_index = 0;
          while (lapd_cyl_3x2v_temp[substring_index] != ' ' && lapd_cyl_3x2v_temp[substring_index] != '\n') {
            substring[substring_index] = lapd_cyl_3x2v_temp[substring_index];
            substring_index += 1;
          }

          char *lapd_cyl_3x2v_temp2 = lapd_cyl_3x2v_output;
          int count = 0;
          while (strstr(lapd_cyl_3x2v_temp2, substring) != NULL) {
            lapd_cyl_3x2v_temp2 = strstr(lapd_cyl_3x2v_temp2, substring);

            count += 1;
            lapd_cyl_3x2v_temp2 += 1;
          }
          if (count == 1) {
            lapd_cyl_3x2v_memoryleakcount[i] += 1;
            lapd_cyl_3x2v_memoryleaks[i] = strcat(lapd_cyl_3x2v_memoryleaks[i], substring);
            lapd_cyl_3x2v_memoryleaks[i] = strcat(lapd_cyl_3x2v_memoryleaks[i], " ");
          }
          
          lapd_cyl_3x2v_temp += 1;
        }
      }

      for (int i = 1; i < lapd_cyl_3x2v_counter + 1; i++) {
        printf("Build number: %d\n", i);
        if (i == 1) {
          printf("Update calls: %d\n", lapd_cyl_3x2v_updatecalls[i]);
          printf("Forward-Euler calls: %d\n", lapd_cyl_3x2v_forwardeuler[i]);
          printf("RK stage-2 failures: %d\n", lapd_cyl_3x2v_rk2failures[i]);
          printf("RK stage-3 failures: %d\n", lapd_cyl_3x2v_rk3failures[i]);
          printf("Species RHS time: %f\n", lapd_cyl_3x2v_speciesrhs[i]);
          printf("Species collision RHS time: %f\n", lapd_cyl_3x2v_speciescollisionsrhs[i]);
          printf("Field RHS time: %f\n", lapd_cyl_3x2v_fieldrhs[i]);
          printf("Species collisional moments time: %f\n", lapd_cyl_3x2v_speciescollisionalmoments[i]);
          printf("Total update time: %f\n", lapd_cyl_3x2v_totalupdate[i]);
          if (lapd_cyl_3x2v_memoryleakcount[i] != 0) {
            printf("Memory leaks: " ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", lapd_cyl_3x2v_memoryleaks[i]);
          }
          else {
            printf("Memory leaks: " ANSI_COLOR_GREEN "None" ANSI_COLOR_RESET "\n");
          }
          printf("Correct: N/A\n\n");
        }
        else {
          if (lapd_cyl_3x2v_updatecalls[i] != lapd_cyl_3x2v_updatecalls[i - 1]) {
            printf("Update calls: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", lapd_cyl_3x2v_updatecalls[i]);
          }
          else {
            printf("Update calls: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", lapd_cyl_3x2v_updatecalls[i]);
          }

          if (lapd_cyl_3x2v_forwardeuler[i] != lapd_cyl_3x2v_forwardeuler[i - 1]) {
            printf("Forward-Euler calls: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", lapd_cyl_3x2v_forwardeuler[i]);
          }
          else {
            printf("Forward-Euler calls: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", lapd_cyl_3x2v_forwardeuler[i]);
          }

          if (lapd_cyl_3x2v_rk2failures[i] > lapd_cyl_3x2v_rk2failures[i - 1]) {
            printf("RK stage-2 failures: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", lapd_cyl_3x2v_rk2failures[i]);
          }
          else {
            printf("RK stage-2 failures: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", lapd_cyl_3x2v_rk2failures[i]);
          }

          if (lapd_cyl_3x2v_rk3failures[i] > lapd_cyl_3x2v_rk3failures[i - 1]) {
            printf("RK stage-3 failures: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", lapd_cyl_3x2v_rk3failures[i]);
          }
          else {
            printf("RK stage-3 failures: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", lapd_cyl_3x2v_rk3failures[i]);
          }

          if (lapd_cyl_3x2v_speciesrhs[i] > lapd_cyl_3x2v_speciesrhs[i - 1]) {
            printf("Species RHS time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", lapd_cyl_3x2v_speciesrhs[i]);
          }
          else {
            printf("Species RHS time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", lapd_cyl_3x2v_speciesrhs[i]);
          }

          if (lapd_cyl_3x2v_speciescollisionsrhs[i] > lapd_cyl_3x2v_speciescollisionsrhs[i - 1]) {
            printf("Species collision RHS time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", lapd_cyl_3x2v_speciescollisionsrhs[i]);
          }
          else {
            printf("Species collision RHS time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", lapd_cyl_3x2v_speciescollisionsrhs[i]);
          }

          if (lapd_cyl_3x2v_fieldrhs[i] > lapd_cyl_3x2v_fieldrhs[i - 1]) {
            printf("Field RHS time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", lapd_cyl_3x2v_fieldrhs[i]);
          }
          else {
            printf("Field RHS time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", lapd_cyl_3x2v_fieldrhs[i]);
          }

          if (lapd_cyl_3x2v_speciescollisionalmoments[i] > lapd_cyl_3x2v_speciescollisionalmoments[i - 1]) {
            printf("Species collisional moments time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", lapd_cyl_3x2v_speciescollisionalmoments[i]);
          }
          else {
            printf("Species collisional moments time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", lapd_cyl_3x2v_speciescollisionalmoments[i]);
          }

          if (lapd_cyl_3x2v_totalupdate[i] > lapd_cyl_3x2v_totalupdate[i - 1]) {
            printf("Total update time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", lapd_cyl_3x2v_totalupdate[i]);
          }
          else {
            printf("Total update time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", lapd_cyl_3x2v_totalupdate[i]);
          }

          if (lapd_cyl_3x2v_memoryleakcount[i] != 0) {
            printf("Memory leaks: " ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", lapd_cyl_3x2v_memoryleaks[i]);
          }
          else {
            printf("Memory leaks: " ANSI_COLOR_GREEN "None" ANSI_COLOR_RESET "\n");
          }

          if ((lapd_cyl_3x2v_updatecalls[i] != lapd_cyl_3x2v_updatecalls[i - 1]) || (lapd_cyl_3x2v_forwardeuler[i] != lapd_cyl_3x2v_forwardeuler[i - 1])
            || (lapd_cyl_3x2v_rk2failures[i] != lapd_cyl_3x2v_rk2failures[i - 1]) || (lapd_cyl_3x2v_rk3failures[i] != lapd_cyl_3x2v_rk3failures[i - 1])) {
            printf("Correct: " ANSI_COLOR_RED "No" ANSI_COLOR_RESET "\n\n");
          }
          else {
            printf("Correct: " ANSI_COLOR_GREEN "Yes" ANSI_COLOR_RESET "\n\n");
          }
        }
      }
    }
  }
  else if (option == 3) {
    system("rm -rf output");
    system("mkdir output");

    FILE *sheath_1x2v_counter_ptr = fopen("output/sheath_1x2v_counter.dat", "w");
    fprintf(sheath_1x2v_counter_ptr, "%d", 1);
    fclose(sheath_1x2v_counter_ptr);

    printf("Running 1x2v Sheath Boundary Test with p = 1...\n");
    system("cd ../; rm -rf ./gk_sheath_1x2v_p1-stat.json");
    system("cd ../; make build/regression/rt_gk_sheath_1x2v_p1 > /dev/null 2>&1");
    char sheath_1x2v_buffer[256];
    snprintf(sheath_1x2v_buffer, 256, "cd ../; ./build/regression/rt_gk_sheath_1x2v_p1 -m > ./ci/output/rt_gk_sheath_1x2v_p1_%d.dat 2>&1", 1);
    system(sheath_1x2v_buffer);
    char sheath_1x2v_buffer2[256];
    snprintf(sheath_1x2v_buffer2, 256, "cd ../; mv ./gk_sheath_1x2v_p1-stat.json ci/output/gk_sheath_1x2v_p1-stat_%d.json", 1);
    system(sheath_1x2v_buffer2);
    printf("Finished 1x2v Sheath Boundary Test with p = 1.\n\n");

    FILE *sheath_2x2v_counter_ptr = fopen("output/sheath_2x2v_counter.dat", "w");
    fprintf(sheath_2x2v_counter_ptr, "%d", 1);
    fclose(sheath_2x2v_counter_ptr);

    printf("Running 2x2v Sheath Boundary Test with p = 1...\n");
    system("cd ../; rm -rf ./gk_sheath_2x2v_p1-stat.json");
    system("cd ../; make build/regression/rt_gk_sheath_2x2v_p1 > /dev/null 2>&1");
    char sheath_2x2v_buffer[256];
    snprintf(sheath_2x2v_buffer, 256, "cd ../; ./build/regression/rt_gk_sheath_2x2v_p1 -m > ./ci/output/rt_gk_sheath_2x2v_p1_%d.dat 2>&1", 1);
    system(sheath_2x2v_buffer);
    char sheath_2x2v_buffer2[256];
    snprintf(sheath_2x2v_buffer2, 256, "cd ../; mv ./gk_sheath_2x2v_p1-stat.json ci/output/gk_sheath_2x2v_p1-stat_%d.json", 1);
    system(sheath_2x2v_buffer2);
    printf("Finished 2x2v Sheath Boundary Test with p = 1.\n\n");

    FILE *sheath_3x2v_counter_ptr = fopen("output/sheath_3x2v_counter.dat", "w");
    fprintf(sheath_3x2v_counter_ptr, "%d", 1);
    fclose(sheath_3x2v_counter_ptr);

    printf("Running 3x2v Sheath Boundary Test with p = 1...\n");
    system("cd ../; rm -rf ./gk_sheath_3x2v_p1-stat.json");
    system("cd ../; make build/regression/rt_gk_sheath_3x2v_p1 > /dev/null 2>&1");
    char sheath_3x2v_buffer[256];
    snprintf(sheath_3x2v_buffer, 256, "cd ../; ./build/regression/rt_gk_sheath_3x2v_p1 -m > ./ci/output/rt_gk_sheath_3x2v_p1_%d.dat 2>&1", 1);
    system(sheath_3x2v_buffer);
    char sheath_3x2v_buffer2[256];
    snprintf(sheath_3x2v_buffer2, 256, "cd ../; mv ./gk_sheath_3x2v_p1-stat.json ci/output/gk_sheath_3x2v_p1-stat_%d.json", 1);
    system(sheath_3x2v_buffer2);
    printf("Finished 3x2v Sheath Boundary Test with p = 1.\n\n");

    FILE *lapd_cart_3x2v_counter_ptr = fopen("output/lapd_cart_3x2v_counter.dat", "w");
    fprintf(lapd_cart_3x2v_counter_ptr, "%d", 1);
    fclose(lapd_cart_3x2v_counter_ptr);

    printf("Running 3x2v LAPD Test (in Cartesian coordinates) with p = 1...\n");
    system("cd ../; rm -rf ./gk_lapd_cart_3x2v_p1-stat.json");
    system("cd ../; make build/regression/rt_gk_lapd_cart_3x2v_p1 > /dev/null 2>&1");
    char lapd_cart_3x2v_buffer[256];
    snprintf(lapd_cart_3x2v_buffer, 256, "cd ../; ./build/regression/rt_gk_lapd_cart_3x2v_p1 -m > ./ci/output/rt_gk_lapd_cart_3x2v_p1_%d.dat 2>&1", 1);
    system(lapd_cart_3x2v_buffer);
    char lapd_cart_3x2v_buffer2[256];
    snprintf(lapd_cart_3x2v_buffer2, 256, "cd ../; mv ./gk_lapd_cart_3x2v_p1-stat.json ci/output/gk_lapd_cart_3x2v_p1-stat_%d.json", 1);
    system(lapd_cart_3x2v_buffer2);
    printf("Finished 3x2v LAPD Test (in Cartesian coordinates) with p = 1.\n\n");

    FILE *lapd_cyl_3x2v_counter_ptr = fopen("output/lapd_cyl_3x2v_counter.dat", "w");
    fprintf(lapd_cyl_3x2v_counter_ptr, "%d", 1);
    fclose(lapd_cyl_3x2v_counter_ptr);

    printf("Running 3x2v LAPD Test (in cylindrical coordinates) with p = 1...\n");
    system("cd ../; rm -rf ./gk_lapd_cyl_3x2v_p1-stat.json");
    system("cd ../; make build/regression/rt_gk_lapd_cyl_3x2v_p1 > /dev/null 2>&1");
    char lapd_cyl_3x2v_buffer[256];
    snprintf(lapd_cyl_3x2v_buffer, 256, "cd ../; ./build/regression/rt_gk_lapd_cyl_3x2v_p1 -m > ./ci/output/rt_gk_lapd_cyl_3x2v_p1_%d.dat 2>&1", 1);
    system(lapd_cyl_3x2v_buffer);
    char lapd_cyl_3x2v_buffer2[256];
    snprintf(lapd_cyl_3x2v_buffer2, 256, "cd ../; mv ./gk_lapd_cyl_3x2v_p1-stat.json ci/output/gk_lapd_cyl_3x2v_p1-stat_%d.json", 1);
    system(lapd_cyl_3x2v_buffer2);
    printf("Finished 3x2v LAPD Test (in cylindrical coordinates) with p = 1.\n\n");
  }
  else if (option == 4) {
    printf("Please select the test you wish to run:\n\n");
    printf("1 - 1x2v Sheath Boundary Test with p = 1\n");
    printf("2 - 2x2v Sheath Boundary Test with p = 1\n");
    printf("3 - 3x3v Sheath Boundary Test with p = 1\n");
    printf("4 - 3x2v LAPD Test (in Cartesian coordinates) with p = 1\n");
    printf("5 - 3x2v LAPD Test (in cylindrical coordinates) with p = 1\n");

    int option2;
    scanf("%d", &option2);
    printf("\n");

    if (option2 == 1) {
      int sheath_1x2v_counter = 0;

      FILE *sheath_1x2v_counter_ptr = fopen("output/sheath_1x2v_counter.dat", "r");
      if (sheath_1x2v_counter_ptr != NULL) {
        fscanf(sheath_1x2v_counter_ptr, "%d", &sheath_1x2v_counter);
      }
      fclose(sheath_1x2v_counter_ptr);

      sheath_1x2v_counter += 1;

      sheath_1x2v_counter_ptr = fopen("output/sheath_1x2v_counter.dat", "w");
      fprintf(sheath_1x2v_counter_ptr, "%d", sheath_1x2v_counter);
      fclose(sheath_1x2v_counter_ptr);

      printf("Running 1x2v Sheath Boundary Test with p = 1...\n");
      system("cd ../; rm -rf ./gk_sheath_1x2v_p1-stat.json");
      system("cd ../; make build/regression/rt_gk_sheath_1x2v_p1 > /dev/null 2>&1");
      char sheath_1x2v_buffer[256];
      snprintf(sheath_1x2v_buffer, 256, "cd ../; ./build/regression/rt_gk_sheath_1x2v_p1 -m > ./ci/output/rt_gk_sheath_1x2v_p1_%d.dat 2>&1", sheath_1x2v_counter);
      system(sheath_1x2v_buffer);
      char sheath_1x2v_buffer2[256];
      snprintf(sheath_1x2v_buffer2, 256, "cd ../; mv ./gk_sheath_1x2v_p1-stat.json ci/output/gk_sheath_1x2v_p1-stat_%d.json", sheath_1x2v_counter);
      system(sheath_1x2v_buffer2);
      printf("Finished 1x2v Sheath Boundary Test with p = 1.\n\n");
    }
    else if (option2 == 2) {
      int sheath_2x2v_counter = 0;

      FILE *sheath_2x2v_counter_ptr = fopen("output/sheath_2x2v_counter.dat", "r");
      if (sheath_2x2v_counter_ptr != NULL) {
        fscanf(sheath_2x2v_counter_ptr, "%d", &sheath_2x2v_counter);
      }
      fclose(sheath_2x2v_counter_ptr);

      sheath_2x2v_counter += 1;

      sheath_2x2v_counter_ptr = fopen("output/sheath_2x2v_counter.dat", "w");
      fprintf(sheath_2x2v_counter_ptr, "%d", sheath_2x2v_counter);
      fclose(sheath_2x2v_counter_ptr);

      printf("Running 2x2v Sheath Boundary Test with p = 1...\n");
      system("cd ../; rm -rf ./gk_sheath_2x2v_p1-stat.json");
      system("cd ../; make build/regression/rt_gk_sheath_2x2v_p1 > /dev/null 2>&1");
      char sheath_2x2v_buffer[256];
      snprintf(sheath_2x2v_buffer, 256, "cd ../; ./build/regression/rt_gk_sheath_2x2v_p1 -m > ./ci/output/rt_gk_sheath_2x2v_p1_%d.dat 2>&1", sheath_2x2v_counter);
      system(sheath_2x2v_buffer);
      char sheath_2x2v_buffer2[256];
      snprintf(sheath_2x2v_buffer2, 256, "cd ../; mv ./gk_sheath_2x2v_p1-stat.json ci/output/gk_sheath_2x2v_p1-stat_%d.json", sheath_2x2v_counter);
      system(sheath_2x2v_buffer2);
      printf("Finished 2x2v Sheath Boundary Test with p = 1.\n\n");
    }
    else if (option2 == 3) {
      int sheath_3x2v_counter = 0;

      FILE *sheath_3x2v_counter_ptr = fopen("output/sheath_3x2v_counter.dat", "r");
      if (sheath_3x2v_counter_ptr != NULL) {
        fscanf(sheath_3x2v_counter_ptr, "%d", &sheath_3x2v_counter);
      }
      fclose(sheath_3x2v_counter_ptr);

      sheath_3x2v_counter += 1;

      sheath_3x2v_counter_ptr = fopen("output/sheath_3x2v_counter.dat", "w");
      fprintf(sheath_3x2v_counter_ptr, "%d", sheath_3x2v_counter);
      fclose(sheath_3x2v_counter_ptr);

      printf("Running 3x2v Sheath Boundary Test with p = 1...\n");
      system("cd ../; rm -rf ./gk_sheath_3x2v_p1-stat.json");
      system("cd ../; make build/regression/rt_gk_sheath_3x2v_p1 > /dev/null 2>&1");
      char sheath_3x2v_buffer[256];
      snprintf(sheath_3x2v_buffer, 256, "cd ../; ./build/regression/rt_gk_sheath_3x2v_p1 -m > ./ci/output/rt_gk_sheath_3x2v_p1_%d.dat 2>&1", sheath_3x2v_counter);
      system(sheath_3x2v_buffer);
      char sheath_3x2v_buffer2[256];
      snprintf(sheath_3x2v_buffer2, 256, "cd ../; mv ./gk_sheath_3x2v_p1-stat.json ci/output/gk_sheath_3x2v_p1-stat_%d.json", sheath_3x2v_counter);
      system(sheath_3x2v_buffer2);
      printf("Finished 3x2v Sheath Boundary Test with p = 1.\n\n");
    }
    else if (option2 == 4) {
      int lapd_cart_3x2v_counter = 0;

      FILE *lapd_cart_3x2v_counter_ptr = fopen("output/lapd_cart_3x2v_counter.dat", "r");
      if (lapd_cart_3x2v_counter_ptr != NULL) {
        fscanf(lapd_cart_3x2v_counter_ptr, "%d", &lapd_cart_3x2v_counter);
      }
      fclose(lapd_cart_3x2v_counter_ptr);

      lapd_cart_3x2v_counter += 1;

      lapd_cart_3x2v_counter_ptr = fopen("output/lapd_cart_3x2v_counter.dat", "w");
      fprintf(lapd_cart_3x2v_counter_ptr, "%d", lapd_cart_3x2v_counter);
      fclose(lapd_cart_3x2v_counter_ptr);

      printf("Running 3x2v LAPD Test (in Cartesian coordinates) with p = 1...\n");
      system("cd ../; rm -rf ./gk_lapd_cart_3x2v_p1-stat.json");
      system("cd ../; make build/regression/rt_gk_lapd_cart_3x2v_p1 > /dev/null 2>&1");
      char lapd_cart_3x2v_buffer[256];
      snprintf(lapd_cart_3x2v_buffer, 256, "cd ../; ./build/regression/rt_gk_lapd_cart_3x2v_p1 -m > ./ci/output/rt_gk_lapd_cart_3x2v_p1_%d.dat 2>&1", lapd_cart_3x2v_counter);
      system(lapd_cart_3x2v_buffer);
      char lapd_cart_3x2v_buffer2[256];
      snprintf(lapd_cart_3x2v_buffer2, 256, "cd ../; mv ./gk_lapd_cart_3x2v_p1-stat.json ci/output/gk_lapd_cart_3x2v_p1-stat_%d.json", lapd_cart_3x2v_counter);
      system(lapd_cart_3x2v_buffer2);
      printf("Finished 3x2v LAPD Test (in Cartesian coordinates) with p = 1.\n\n");
    }
    else if (option2 == 5) {
      int lapd_cyl_3x2v_counter = 0;

      FILE *lapd_cyl_3x2v_counter_ptr = fopen("output/lapd_cyl_3x2v_counter.dat", "r");
      if (lapd_cyl_3x2v_counter_ptr != NULL) {
        fscanf(lapd_cyl_3x2v_counter_ptr, "%d", &lapd_cyl_3x2v_counter);
      }
      fclose(lapd_cyl_3x2v_counter_ptr);

      lapd_cyl_3x2v_counter += 1;

      lapd_cyl_3x2v_counter_ptr = fopen("output/lapd_cyl_3x2v_counter.dat", "w");
      fprintf(lapd_cyl_3x2v_counter_ptr, "%d", lapd_cyl_3x2v_counter);
      fclose(lapd_cyl_3x2v_counter_ptr);

      printf("Running 3x2v LAPD Test (in cylindrical coordinates) with p = 1...\n");
      system("cd ../; rm -rf ./gk_lapd_cyl_3x2v_p1-stat.json");
      system("cd ../; make build/regression/rt_gk_lapd_cyl_3x2v_p1 > /dev/null 2>&1");
      char lapd_cyl_3x2v_buffer[256];
      snprintf(lapd_cyl_3x2v_buffer, 256, "cd ../; ./build/regression/rt_gk_lapd_cyl_3x2v_p1 -m > ./ci/output/rt_gk_lapd_cyl_3x2v_p1_%d.dat 2>&1", lapd_cyl_3x2v_counter);
      system(lapd_cyl_3x2v_buffer);
      char lapd_cyl_3x2v_buffer2[256];
      snprintf(lapd_cyl_3x2v_buffer2, 256, "cd ../; mv ./gk_lapd_cyl_3x2v_p1-stat.json ci/output/gk_lapd_cyl_3x2v_p1-stat_%d.json", lapd_cyl_3x2v_counter);
      system(lapd_cyl_3x2v_buffer2);
      printf("Finished 3x2v LAPD Test (in cylindrical coordinates) with p = 1.\n\n");
    }
  }

  return 0;
}