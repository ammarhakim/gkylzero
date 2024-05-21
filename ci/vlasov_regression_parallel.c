#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_YELLOW "\x1b[33m"
#define ANSI_COLOR_BLUE "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN "\x1b[36m"
#define ANSI_COLOR_RESET "\x1b[0m"

#define RELATIVE_TOLERANCE pow(10.0, -16.0)

int system(const char *command);

void
runTestParallel(const char* test_name, const char* test_name_human, const int test_dimensions, const int test_cuts,
  const int test_output_count, const char test_outputs[][64])
{
  int counter = 0;

  char counter_buffer[128];
  snprintf(counter_buffer, 128, "ci/output_parallel/%s_counter.dat", test_name);
  FILE *counter_ptr = fopen(counter_buffer, "r");
  if (counter_ptr != NULL) {
    fscanf(counter_ptr, "%d", &counter);
    fclose(counter_ptr);
  }

  counter += 1;

  counter_ptr = fopen(counter_buffer, "w");
  fprintf(counter_ptr, "%d", counter);
  fclose(counter_ptr);

  printf("Running %s...\n", test_name_human);

  char command_buffer1[256];
  snprintf(command_buffer1, 256, "rm -rf ./%s-stat.json", test_name);
  system(command_buffer1);
  
  char command_buffer2[256];
  snprintf(command_buffer2, 256, "make cuda-build/regression/rt_%s > /dev/null 2>&1", test_name);
  system(command_buffer2);

  char command_buffer3[256];
  if (test_dimensions == 1) {
    snprintf(command_buffer3, 256, "mpirun -np %d ./cuda-build/regression/rt_%s -m -M -c %d > ./ci/output_parallel/rt_%s_%d.dat 2>&1",
      test_cuts, test_name, test_cuts, test_name, counter);
  }
  else if (test_dimensions == 2) {
    snprintf(command_buffer3, 256, "mpirun -np %d ./cuda-build/regression/rt_%s -m -M -d %d > ./ci/output_parallel/rt_%s_%d.dat 2>&1",
      test_cuts, test_name, test_cuts, test_name, counter);
  }
  else if (test_dimensions == 3) {
    snprintf(command_buffer3, 256, "mpirun -np %d ./cuda-build/regression/rt_%s -m -M -e %d > ./ci/output_parallel/rt_%s_%d.dat 2>&1",
      test_cuts, test_name, test_cuts, test_name, counter);
  }
  system(command_buffer3);

  char file_buffer1[128];
  snprintf(file_buffer1, 128, "./%s-stat.json", test_name);
  FILE *file_ptr1 = fopen(file_buffer1, "r");
  if (file_ptr1 == NULL) {
    printf("*** Something catastrophic happened. Test aborting... ***\n");
  }
  else {
    char command_buffer4[256];
    snprintf(command_buffer4, 256, "mv ./%s-stat.json ci/output_parallel/%s-stat_%d.json", test_name, test_name, counter);
    system(command_buffer4);
  }

  for (int i = 0; i < test_output_count; i++) {
    char file_buffer2[128];
    snprintf(file_buffer2, 128, "./%s-%s.gkyl", test_name, test_outputs[i]);
    FILE *file_ptr2 = fopen(file_buffer2, "r");
    if (file_ptr2 == NULL) {
      printf("*** Something catastrophic happened. Test aborting... ***\n");
    }
    else {
      char command_buffer5[256];
      snprintf(command_buffer5, 256, "mv ./%s-%s.gkyl ci/output_parallel/%s-%s_%d.gkyl", test_name, test_outputs[i], test_name, test_outputs[i], counter);
      system(command_buffer5);
    }
  }

  printf("Finished %s.\n\n", test_name_human);
}

void
analyzeTestOutputParallel(const char* test_name, const char* test_name_human, const int test_output_count, const char test_outputs[][64])
{
  printf("%s:\n\n", test_name_human);

  int counter = 0;

  char counter_buffer[64];
  snprintf(counter_buffer, 64, "ci/output_parallel/%s_counter.dat", test_name);
  FILE *counter_ptr = fopen(counter_buffer, "r");
  if (counter_ptr != NULL) {
    fscanf(counter_ptr, "%d", &counter);
    fclose(counter_ptr);
  }

  int failure[counter + 1];
  int updatecalls[counter + 1];
  int forwardeuler[counter + 1];
  int rk2failures[counter + 1];
  int rk3failures[counter + 1];
  double speciesrhs[counter + 1];
  double speciescollisionsrhs[counter + 1];
  double fieldrhs[counter + 1];
  double speciescollisionalmoments[counter + 1];
  double totalupdate[counter + 1];
  int memoryleakcount[counter + 1];
  char *memoryleaks[counter + 1];
  long double averages[counter + 1][test_output_count];

  for (int i = 1; i < counter + 1; i++) {
    char *output;
    long file_size;
    char buffer[128];
    snprintf(buffer, 128, "ci/output_parallel/rt_%s_%d.dat", test_name, i);

    FILE *output_ptr = fopen(buffer, "rb");
    fseek(output_ptr, 0, SEEK_END);
    file_size = ftell(output_ptr);
    rewind(output_ptr);
    output = calloc(file_size, (sizeof(char)));
    fread(output, sizeof(char), file_size, output_ptr);
    fclose(output_ptr);

    failure[i] = 0;

    updatecalls[i] = 0;
    if (strstr(output, "Number of update calls ") != NULL) {
      char *full_substring = strstr(output, "Number of update calls ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Number of update calls ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Number of update calls ")];
        substring_index += 1;
      }

      char *end_ptr;
      updatecalls[i] = strtol(substring, &end_ptr, 10);
    }
    else {
      failure[i] = 1;
    }

    forwardeuler[i] = 0;
    if (strstr(output, "Number of forward-Euler calls ") != NULL) {
      char *full_substring = strstr(output, "Number of forward-Euler calls ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Number of forward-Euler calls ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Number of forward-Euler calls ")];
        substring_index += 1;
      }

      char *end_ptr;
      forwardeuler[i] = strtol(substring, &end_ptr, 10);
    }
    else {
      failure[i] = 1;
    }

    rk2failures[i] = 0;
    if (strstr(output, "Number of RK stage-2 failures ") != NULL) {
      char *full_substring = strstr(output, "Number of RK stage-2 failures ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Number of RK stage-2 failures ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Number of RK stage-2 failures ")];
        substring_index += 1;
      }

      char *end_ptr;
      rk2failures[i] = strtol(substring, &end_ptr, 10);
    }
    else {
      failure[i] = 1;
    }

    rk3failures[i] = 0;
    if (strstr(output, "Number of RK stage-3 failures ") != NULL) {
      char *full_substring = strstr(output, "Number of RK stage-3 failures ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Number of RK stage-3 failures ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Number of RK stage-3 failures ")];
        substring_index += 1;
      }

      char *end_ptr;
      rk3failures[i] = strtol(substring, &end_ptr, 10);
    }
    else {
      failure[i] = 1;
    }

    speciesrhs[i] = 0.0;
    if (strstr(output, "Species RHS calc took ") != NULL) {
      char *full_substring = strstr(output, "Species RHS calc took ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Species RHS calc took ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Species RHS calc took ")];
        substring_index += 1;
      }

      char *end_ptr;
      speciesrhs[i] = strtod(substring, &end_ptr);
    }
    else {
      failure[i] = 1;
    }

    speciescollisionsrhs[i] = 0.0;
    if (strstr(output, "Species collisions RHS calc took ") != NULL) {
      char *full_substring = strstr(output, "Species collisions RHS calc took ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Species collisions RHS calc took ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Species collisions RHS calc took ")];
        substring_index += 1;
      }

      char *end_ptr;
      speciescollisionsrhs[i] = strtod(substring, &end_ptr);
    }
    else {
      failure[i] = 1;
    }

    fieldrhs[i] = 0.0;
    if (strstr(output, "Field RHS calc took ") != NULL) {
      char *full_substring = strstr(output, "Field RHS calc took ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Field RHS calc took ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Field RHS calc took ")];
        substring_index += 1;
      }

      char *end_ptr;
      fieldrhs[i] = strtod(substring, &end_ptr);
    }
    else {
      failure[i] = 1;
    }

    speciescollisionalmoments[i] = 0.0;
    if (strstr(output, "Species collisional moments took ") != NULL) {
      char *full_substring = strstr(output, "Species collisional moments took ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Species collisional moments took ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Species collisional moments took ")];
        substring_index += 1;
      }

      char *end_ptr;
      speciescollisionalmoments[i] = strtod(substring, &end_ptr);
    }
    else {
      failure[i] = 1;
    }

    totalupdate[i] = 0.0;
    if (strstr(output, "Total updates took ") != NULL) {
      char *full_substring = strstr(output, "Total updates took ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Total updates took ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Total updates took ")];
        substring_index += 1;
      }

      char *end_ptr;
      totalupdate[i] = strtod(substring, &end_ptr);
    }
    else {
      failure[i] = 1;
    }

    if (failure[i] == 0) {
      char *temp = output;
      memoryleakcount[i] = 0;
      memoryleaks[i] = (char*)calloc(8192, sizeof(char));
      while (strstr(temp, "0x") != NULL) {
        temp = strstr(temp, "0x");

        char substring[64];
        for (int j = 0; j < 64; j++) {
          substring[j] = '\0';
        }

        int substring_index = 0;
        int valid_substring = 1;
        while (temp[substring_index] != ' ' && temp[substring_index] != '\n') {
          if (temp[substring_index] != '0' && temp[substring_index] != '1' && temp[substring_index] != '2' && temp[substring_index] != '3' && temp[substring_index] != '4'
            && temp[substring_index] != '5' && temp[substring_index] != '6' && temp[substring_index] != '7' && temp[substring_index] != '8' && temp[substring_index] != '9'
            && temp[substring_index] != 'a' && temp[substring_index] != 'b' && temp[substring_index] != 'c' && temp[substring_index] != 'd' && temp[substring_index] != 'e'
            && temp[substring_index] != 'f' && temp[substring_index] != 'x') {
            valid_substring = 0;
          }

          substring[substring_index] = temp[substring_index];
          substring_index += 1;
        }

        char *temp2 = output;
        int count = 0;
        while (strstr(temp2, substring) != NULL) {
          temp2 = strstr(temp2, substring);

          count += 1;
          temp2 += 1;
        }
        if (count == 1 && valid_substring == 1) {
          memoryleakcount[i] += 1;
          memoryleaks[i] = strcat(memoryleaks[i], substring);
          memoryleaks[i] = strcat(memoryleaks[i], " ");
        }
        
        temp += 1;
      }
    }

    for (int j = 0; j < test_output_count; j++) {
      char *data;
      long data_file_size;
      char data_buffer[256];
      snprintf(data_buffer, 256, "ci/output_parallel/%s-%s_%d.gkyl", test_name, test_outputs[j], i);

      FILE *data_ptr = fopen(data_buffer, "rb");

      if (data_ptr == NULL) {
        failure[i] = 1;
      }
      else {
        fseek(data_ptr, 0, SEEK_END);
        data_file_size = ftell(data_ptr);
        rewind(data_ptr);
        data = calloc(data_file_size, (sizeof(char)));
        fread(data, sizeof(char), data_file_size, data_ptr);
        fclose(data_ptr);

        long long total = 0;
        for (long k = 0; k < data_file_size; k++) {
          total += (long long)abs((int)data[k]);
        }

        averages[i][j] = (long double)total / (long double)data_file_size;
      }
    }
  }

  for (int i = 1; i < counter + 1; i++) {
    printf("Build number: %d\n", i);

    if (failure[i] == 1) {
      printf(ANSI_COLOR_RED "*** Catastrophic test failure ***" ANSI_COLOR_RESET "\n\n");
    }
    else {
      if (i == 1 || failure[i - 1] == 1) {
        printf("Update calls: %d\n", updatecalls[i]);
        printf("Forward-Euler calls: %d\n", forwardeuler[i]);
        printf("RK stage-2 failures: %d\n", rk2failures[i]);
        printf("RK stage-3 failures: %d\n", rk3failures[i]);
        printf("Species RHS time: %f\n", speciesrhs[i]);
        printf("Species collision RHS time: %f\n", speciescollisionsrhs[i]);
        printf("Field RHS time: %f\n", fieldrhs[i]);
        printf("Species collisional moments time: %f\n", speciescollisionalmoments[i]);
        printf("Total update time: %f\n", totalupdate[i]);
        if (memoryleakcount[i] != 0) {
          printf("Memory leaks: " ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", memoryleaks[i]);
        }
        else {
          printf("Memory leaks: " ANSI_COLOR_GREEN "None" ANSI_COLOR_RESET "\n");
        }
        printf("Correct: N/A\n\n");
      }
      else {
        if (updatecalls[i] != updatecalls[i - 1]) {
          printf("Update calls: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", updatecalls[i]);
        }
        else {
          printf("Update calls: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", updatecalls[i]);
        }

        if (forwardeuler[i] != forwardeuler[i - 1]) {
          printf("Forward-Euler calls: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", forwardeuler[i]);
        }
        else {
          printf("Forward-Euler calls: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", forwardeuler[i]);
        }

        if (rk2failures[i] > rk2failures[i - 1]) {
          printf("RK stage-2 failures: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", rk2failures[i]);
        }
        else {
          printf("RK stage-2 failures: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", rk2failures[i]);
        }

        if (rk3failures[i] > rk3failures[i - 1]) {
          printf("RK stage-3 failures: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", rk3failures[i]);
        }
        else {
          printf("RK stage-3 failures: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", rk3failures[i]);
        }

        if (speciesrhs[i] > speciesrhs[i - 1]) {
          if (speciesrhs[i - 1] > pow(10.0, -8.0)) {
            printf("Species RHS time: " ANSI_COLOR_RED "%f (+%.2f%%)" ANSI_COLOR_RESET "\n", speciesrhs[i],
              (((double)speciesrhs[i] / (double)speciesrhs[i - 1]) - 1.0) * 100.0);
          } else {
            printf("Species RHS time: " ANSI_COLOR_RED "%f (N/A)" ANSI_COLOR_RESET "\n", speciesrhs[i]);
          }
        }
        else {
          if (speciesrhs[i - 1] > pow(10.0, -8.0)) {
            printf("Species RHS time: " ANSI_COLOR_GREEN "%f (%.2f%%)" ANSI_COLOR_RESET "\n", speciesrhs[i],
              (((double)speciesrhs[i] / (double)speciesrhs[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Species RHS time: " ANSI_COLOR_GREEN "%f (N/A)" ANSI_COLOR_RESET "\n", speciesrhs[i]);
          }
        }

        if (speciescollisionsrhs[i] > speciescollisionsrhs[i - 1]) {
          if (speciescollisionsrhs[i - 1] > pow(10.0, -8.0)) {
            printf("Species collision RHS time: " ANSI_COLOR_RED "%f (+%.2f%%)" ANSI_COLOR_RESET "\n", speciescollisionsrhs[i],
              (((double)speciescollisionsrhs[i] / (double)speciescollisionsrhs[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Species collision RHS time: " ANSI_COLOR_RED "%f (N/A)" ANSI_COLOR_RESET, speciescollisionsrhs[i]);
          }
        }
        else {
          if (speciescollisionsrhs[i - 1] > pow(10.0, -8.0)) {
            printf("Species collision RHS time: " ANSI_COLOR_GREEN "%f (%.2f%%)" ANSI_COLOR_RESET "\n", speciescollisionsrhs[i],
              (((double)speciescollisionsrhs[i] / (double)speciescollisionsrhs[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Species collision RHS time: " ANSI_COLOR_GREEN "%f (N/A)" ANSI_COLOR_RESET "\n", speciescollisionsrhs[i]);
          }
        }

        if (fieldrhs[i] > fieldrhs[i - 1]) {
          if (fieldrhs[i - 1] > pow(10.0, -8.0)) {
            printf("Field RHS time: " ANSI_COLOR_RED "%f (+%.2f%%)" ANSI_COLOR_RESET "\n", fieldrhs[i],
              (((double)fieldrhs[i] / (double)fieldrhs[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Field RHS time: " ANSI_COLOR_RED "%f (N/A)" ANSI_COLOR_RESET "\n", fieldrhs[i]);
          }
        }
        else {
          if (fieldrhs[i - 1] > pow(10.0, -8.0)) {
            printf("Field RHS time: " ANSI_COLOR_GREEN "%f (%.2f%%)" ANSI_COLOR_RESET "\n", fieldrhs[i],
              (((double)fieldrhs[i] / (double)fieldrhs[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Field RHS time: " ANSI_COLOR_GREEN "%f (N/A)" ANSI_COLOR_RESET "\n", fieldrhs[i]);
          }
        }

        if (speciescollisionalmoments[i] > speciescollisionalmoments[i - 1]) {
          if (speciescollisionalmoments[i - 1] > pow(10.0, -8.0)) {
            printf("Species collisional moments time: " ANSI_COLOR_RED "%f (+%.2f%%)" ANSI_COLOR_RESET "\n", speciescollisionalmoments[i],
              (((double)speciescollisionalmoments[i] / (double)speciescollisionalmoments[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Species collisional moments time: " ANSI_COLOR_RED "%f (N/A)" ANSI_COLOR_RESET "\n", speciescollisionalmoments[i]);
          }
        }
        else {
          if (speciescollisionalmoments[i - 1] > pow(10.0, -8.0)) {
            printf("Species collisional moments time: " ANSI_COLOR_GREEN "%f (%.2f%%)" ANSI_COLOR_RESET "\n", speciescollisionalmoments[i],
              (((double)speciescollisionalmoments[i] / (double)speciescollisionalmoments[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Species collisional moments time: " ANSI_COLOR_GREEN "%f (N/A)" ANSI_COLOR_RESET "\n", speciescollisionalmoments[i]);
          }
        }

        if (totalupdate[i] > totalupdate[i - 1]) {
          if (totalupdate[i  - 1] > pow(10.0, -8.0)) {
            printf("Total update time: " ANSI_COLOR_RED "%f (+%.2f%%)" ANSI_COLOR_RESET "\n", totalupdate[i],
              (((double)totalupdate[i] / (double)totalupdate[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Total update time: " ANSI_COLOR_RED "%f (N/A)" ANSI_COLOR_RESET "\n", totalupdate[i]);
          }
        }
        else {
          if (totalupdate[i - 1] > pow(10.0, -8.0)) {
            printf("Total update time: " ANSI_COLOR_GREEN "%f (%.2f%%)" ANSI_COLOR_RESET "\n", totalupdate[i],
              (((double)totalupdate[i] / (double)totalupdate[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Total update time: " ANSI_COLOR_GREEN "%f (N/A)" ANSI_COLOR_RESET "\n", totalupdate[i]);
          }
        }

        if (memoryleakcount[i] != 0) {
          printf("Memory leaks: " ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", memoryleaks[i]);
        }
        else {
          printf("Memory leaks: " ANSI_COLOR_GREEN "None" ANSI_COLOR_RESET "\n");
        }

        int correct = 1;
        for (int j = 0; j < test_output_count; j++) {
          if (fabsl(averages[i][j] - averages[i - 1][j]) > RELATIVE_TOLERANCE) {
            correct = 0;
          }
        }

        if ((updatecalls[i] != updatecalls[i - 1]) || (forwardeuler[i] != forwardeuler[i - 1]) ||
          (rk2failures[i] != rk2failures[i - 1]) || (rk3failures[i] != rk3failures[i - 1]) || (correct != 1)) {
          printf("Correct: " ANSI_COLOR_RED "No" ANSI_COLOR_RESET "\n\n");
        }
        else {
          printf("Correct: " ANSI_COLOR_GREEN "Yes" ANSI_COLOR_RESET "\n\n");
        }
      }
    }
  }
}

void
regenerateTestParallel(const char* test_name, const int test_output_count, const char test_outputs[][64])
{
  int counter = 0;

  char counter_buffer[64];
  snprintf(counter_buffer, 64, "ci/output_parallel/%s_counter.dat", test_name);
  FILE *counter_ptr = fopen(counter_buffer, "r");
  if (counter_ptr != NULL) {
    fscanf(counter_ptr, "%d", &counter);
  }
  fclose(counter_ptr);

  for (int i = 1 ; i < counter + 1; i++) {
    char command_buffer[128];
    snprintf(command_buffer, 128, "rm -rf ci/output_parallel/rt_%s_%d.dat", test_name, i);
    system(command_buffer);

    char command_buffer2[128];
    snprintf(command_buffer2, 128, "rm -rf ci/output_parallel/%s-stat_%d.json", test_name, i);
    system(command_buffer2);

    for (int j = 0; j < test_output_count; j++) {
      char command_buffer3[256];
      snprintf(command_buffer3, 256, "rm -rf ci/output_parallel/%s-%s_%d.gkyl", test_name, test_outputs[j], i);
      system(command_buffer3);
    }
  }

  counter_ptr = fopen(counter_buffer, "w");
  fprintf(counter_ptr, "%d", 0);
  fclose(counter_ptr);
}

int
main(int argc, char **argv)
{
  int test_count = 32;
  char test_names[32][64] = {
    "dg_accel_1x1v",
    "dg_euler_sodshock_p1",
    "dg_euler_sodshock_p2",
    "dg_euler_p_perturbation_p1",
    "dg_euler_p_perturbation_p2",
    "dg_euler_kh_2d",
    "vlasov_lbo_cross_1x1v_p2",
    "vlasov_lbo_cross_1x2v_p2",
    "vlasov_lbo_relax_1x1v_p1",
    "vlasov_lbo_relax_1x1v_p2",
    "vlasov_lbo_relax_1x2v_p1",
    "vlasov_lbo_relax_1x2v_p2",
    "vlasov_lbo_relax_1x3v_p1",
    "vlasov_lbo_relax_1x3v_p2",
    "vlasov_bgk_relax_1x1v_p1",
    "vlasov_bgk_relax_1x1v_p2",
    "vlasov_bgk_relax_1x2v_p1",
    "vlasov_bgk_relax_1x2v_p2",
    "vlasov_bgk_relax_1x3v_p1",
    "vlasov_bgk_relax_1x3v_p2",
    "vlasov_sr_bgk_relax_1x1v_p2",
    "vlasov_neut_bgk_sodshock_1x1v_p1",
    "vlasov_neut_bgk_sodshock_1x1v_p2",
    "vlasov_neut_bgk_sodshock_1x2v_p1",
    "vlasov_neut_bgk_sodshock_1x2v_p2",
    "vlasov_neut_bgk_sodshock_1x3v_p1",
    "vlasov_neut_bgk_sodshock_1x3v_p2",
    "vlasov_sr_neut_bgk_sodshock_1x1v_p2",
    "vlasov_neut_lbo_sodshock_1x1v_p2",
    "vlasov_neut_lbo_sodshock_1x2v_p2",
    "vlasov_neut_lbo_sodshock_1x3v_p2",
    "vlasov_neut_lbo_wall",
  };
  char test_names_human[32][128] = {
    "1x1v Acceleration Test with p = 1",
    "Euler Sod-Type Shock Tube Test p = 1",
    "Euler Sod-Type Shock Tube Test p = 2",
    "Euler Pressure Perturbation Test with p = 1",
    "Euler Pressure Perturbation Test with p = 2",
    "Euler 2D Kelvin-Helmholtz Instability Test with p = 1",
    "1x1v LBO Cross Collision Neutrals Test with p = 2",
    "1x2v LBO Cross Collision Neutrals Test with p = 2",
    "1x1v LBO Collision Relaxation Test with p = 1",
    "1x1v LBO Collision Relaxation Test with p = 2",
    "1x2v LBO Collision Relaxation Test with p = 1",
    "1x2v LBO Collision Relaxation Test with p = 2",
    "1x3v LBO Collision Relaxation Test with p = 1",
    "1x3v LBO Collision Relaxation Test with p = 2",
    "1x1v BGK Collision Relaxation Test with p = 1",
    "1x1v BGK Collision Relaxation Test with p = 2",
    "1x2v BGK Collision Relaxation Test with p = 1",
    "1x2v BGK Collision Relaxation Test with p = 2",
    "1x3v BGK Collision Relaxation Test with p = 1",
    "1x3v BGK Collision Relaxation Test with p = 2",
    "1x1v BGK Relativistic Collision Relaxation Test with p = 1",
    "1x1v BGK Collision Neutrals Shock Tube Test with p = 1",
    "1x1v BGK Collision Neutrals Shock Tube Test with p = 2",
    "1x2v BGK Collision Neutrals Shock Tube Test with p = 1",
    "1x2v BGK Collision Neutrals Shock Tube Test with p = 2",
    "1x3v BGK Collision Neutrals Shock Tube Test with p = 1",
    "1x3v BGK Collision Neutrals Shock Tube Test with p = 2",
    "1x1v BGK Relativistic Collision Neutrals Shock Tube Test with p = 2",
    "1x1v LBO Collision Neutrals Shock Tube Test with p = 2",
    "1x2v LBO Collision Neutrals Shock Tube Test with p = 2",
    "1x3v LBO Collision Neutrals Shock Tube Test with p = 2",
    "1x1v LBO Collision Wall Boundary Test with p = 2",
  };
  int test_dimensions[32] = { 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
  int test_cuts[32] = { 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 };
  int test_output_count[32] = { 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
  char test_outputs[32][64][64] = {
    { "elc_1" },
    { "euler_1" },
    { "euler_1" },
    { "euler_1" },
    { "euler_1" },
    { "euler_1" },
    { "neut1_1", "neut2_1" },
    { "neut1_1", "neut2_1" },
    { "bump_1", "square_1" },
    { "bump_1", "square_1" },
    { "bump_1", "square_1" },
    { "bump_1", "square_1" },
    { "bump_1", "square_1" },
    { "bump_1", "square_1" },
    { "bump_1", "square_1" },
    { "bump_1", "square_1" },
    { "bump_1", "square_1" },
    { "bump_1", "square_1" },
    { "bump_1", "square_1" },
    { "bump_1", "square_1" },
    { "bump_1", "square_1" },
    { "neut_1" },
    { "neut_1" },
    { "neut_1" },
    { "neut_1" },
    { "neut_1" },
    { "neut_1" },
    { "neut_1" },
    { "neut_1" },
    { "neut_1" },
    { "neut_1" },
    { "neut_1" },
  };

  system("clear");
  system("mkdir -p ci/output_parallel");

  printf("** Gkeyll Vlasov Automated Regression System (Parallel Version) **\n\n");

  if (argc > 1) {
    char *arg_ptr;

    if (strtol(argv[1], &arg_ptr, 10) == 1) {
      for (int i = 0; i < test_count; i++) {
        runTestParallel(test_names[i], test_names_human[i], test_dimensions[i], test_cuts[i], test_output_count[i], test_outputs[i]);
      }
    }
    else if (strtol(argv[1], &arg_ptr, 10) == 2) {
      for (int i = 0; i < test_count; i++) {
        analyzeTestOutputParallel(test_names[i], test_names_human[i], test_output_count[i], test_outputs[i]);
      }
    }
    else if (strtol(argv[1], &arg_ptr, 10) == 3) {
      if (argc > 2) {
        if (strtol(argv[2], &arg_ptr, 10) >= 1 && strtol(argv[2], &arg_ptr, 10) <= test_count) {
          runTestParallel(test_names[strtol(argv[2], &arg_ptr, 10) - 1], test_names_human[strtol(argv[2], &arg_ptr, 10) - 1],
            test_dimensions[strtol(argv[2], &arg_ptr, 10) - 1], test_cuts[strtol(argv[2], &arg_ptr, 10) - 1],
            test_output_count[strtol(argv[2], &arg_ptr, 10) - 1], test_outputs[strtol(argv[2], &arg_ptr, 10) - 1]);
        }
        else {
          printf("Invalid test!\n");
        }
      }
      else {
        printf("Must specify which test to run!\n");
      }
    }
    else if (strtol(argv[1], &arg_ptr, 10) == 4) {
      if (argc > 2) {
        if (strtol(argv[2], &arg_ptr, 10) >= 1 && strtol(argv[2], &arg_ptr, 10) <= test_count) {
          analyzeTestOutputParallel(test_names[strtol(argv[2], &arg_ptr, 10) - 1], test_names_human[strtol(argv[2], &arg_ptr, 10) - 1],
            test_output_count[strtol(argv[2], &arg_ptr, 10) - 1], test_outputs[strtol(argv[2], &arg_ptr, 10) - 1]);
        }
        else {
          printf("Invalid test!\n");
        }
      }
      else {
        printf("Must specify which test results to view!\n");
      }
    }
    else if (strtol(argv[1], &arg_ptr, 10) == 5) {
      for (int i = 0; i < test_count; i++) {
        regenerateTestParallel(test_names[i], test_output_count[i], test_outputs[i]);
        runTestParallel(test_names[i], test_names_human[i], test_dimensions[i], test_cuts[i], test_output_count[i], test_outputs[i]);
      }
    }
    else if (strtol(argv[1], &arg_ptr, 10) == 6) {
      if (argc > 2) {
        if (strtol(argv[2], &arg_ptr, 10) >= 1 && strtol(argv[2], &arg_ptr, 10) <= test_count) {
          regenerateTestParallel(test_names[strtol(argv[2], &arg_ptr, 10) - 1], test_output_count[strtol(argv[2], &arg_ptr, 10) - 1],
            test_outputs[strtol(argv[2], &arg_ptr, 10) - 1]);
          runTestParallel(test_names[strtol(argv[2], &arg_ptr, 10) - 1], test_names_human[strtol(argv[2], &arg_ptr, 10) - 1],
            test_dimensions[strtol(argv[2], &arg_ptr, 10) - 1], test_cuts[strtol(argv[2], &arg_ptr, 10) - 1],
            test_output_count[strtol(argv[2], &arg_ptr, 10) - 1], test_outputs[strtol(argv[2], &arg_ptr, 10) - 1]);
        }
        else {
          printf("Invalid test!\n");
        }
      }
      else {
        printf("Must specify which test results to (re)generate!\n");
      }
    }
  }
  else {
    while (1) {
      printf("Please select an option to proceed:\n\n");
      printf("1 - Run Full Regression Suite\n");
      printf("2 - View All Regression Results\n");
      printf("3 - Run Specific Regression Test\n");
      printf("4 - View Specific Regression Result\n");
      printf("5 - (Re)generate All Accepted Results\n");
      printf("6 - (Re)generate Specific Accepted Result\n");
      printf("7 - Exit\n");

      int option;
      scanf("%d", &option);
      printf("\n");

      if (option == 1) {
        for (int i = 0; i < test_count; i++) {
          runTestParallel(test_names[i], test_names_human[i], test_dimensions[i], test_cuts[i], test_output_count[i], test_outputs[i]);
        }
      }
      else if (option == 2) {
        for (int i = 0; i < test_count; i++) {
          analyzeTestOutputParallel(test_names[i], test_names_human[i], test_output_count[i], test_outputs[i]);
        }
      }
      else if (option == 3) {
        printf("Please select the test you wish to run:\n\n");
        for (int i = 0; i < test_count; i++) {
          printf("%d - %s\n", i + 1, test_names_human[i]);
        }

        int option2;
        scanf("%d", &option2);
        printf("\n");

        if (option2 >= 1 && option2 <= test_count) {
          runTestParallel(test_names[option2 - 1], test_names_human[option2 - 1], test_dimensions[option2 - 1], test_cuts[option2 - 1],
            test_output_count[option2 - 1], test_outputs[option2 - 1]);
        }
        else {
          printf("Invalid test!\n\n");
        }
      }
      else if (option == 4) {
        printf("Please select the test whose results you wish to view:\n\n");
        for (int i = 0; i < test_count; i++) {
          printf("%d - %s\n", i + 1, test_names_human[i]);
        }

        int option2;
        scanf("%d", &option2);
        printf("\n");

        if (option2 >= 1 && option2 <= test_count) {
          analyzeTestOutputParallel(test_names[option2 - 1], test_names_human[option2 - 1], test_output_count[option2 - 1], test_outputs[option2 - 1]);
        }
        else {
          printf("Invalid test!\n\n");
        }
      }
      else if (option == 5) {
        for (int i = 0; i < test_count; i++) {
          regenerateTestParallel(test_names[i], test_output_count[i], test_outputs[i]);
          runTestParallel(test_names[i], test_names_human[i], test_dimensions[i], test_cuts[i], test_output_count[i], test_outputs[i]);
        }
      }
      else if (option == 6) {
        printf("Please select the test whose accepted result you wish to (re)generate:\n\n");
        for (int i = 0; i < test_count; i++) {
          printf("%d - %s\n", i + 1, test_names_human[i]);
        }

        int option2;
        scanf("%d", &option2);
        printf("\n");

        if (option2 >= 1 && option2 <= test_count) {
          regenerateTestParallel(test_names[option2 - 1], test_output_count[option2 - 1], test_outputs[option2 - 1]);
          runTestParallel(test_names[option2 - 1], test_names_human[option2 - 1], test_dimensions[option2 - 1], test_cuts[option2 - 1],
            test_output_count[option2 - 1], test_outputs[option2 - 1]);
        }
        else {
          printf("Invalid test!\n\n");
        }
      }
      else if (option == 7) {
        break;
      }
      else {
        printf("Invalid selection!\n\n");
      }
    }
  }

  return 0;
}