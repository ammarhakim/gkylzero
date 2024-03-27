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

void
runTest(const char* test_name, const char* test_name_human)
{
  int counter = 0;

  char counter_buffer[64];
  snprintf(counter_buffer, 64, "output/%s_counter.dat", test_name);
  FILE *counter_ptr = fopen(counter_buffer, "r");
  if (counter_ptr != NULL) {
    fscanf(counter_ptr, "%d", &counter);
  }
  fclose(counter_ptr);

  counter += 1;

  counter_ptr = fopen(counter_buffer, "w");
  fprintf(counter_ptr, "%d", counter);
  fclose(counter_ptr);

  printf("Running %s...\n", test_name_human);

  char command_buffer1[128];
  snprintf(command_buffer1, 128, "cd ../; rm -rf ./%s-stat.json", test_name);
  system(command_buffer1);
  
  char command_buffer2[128];
  snprintf(command_buffer2, 128, "cd ../; make build/regression/rt_%s > /dev/null 2>&1", test_name);
  system(command_buffer2);

  char command_buffer3[256];
  snprintf(command_buffer3, 256, "cd ../; ./build/regression/rt_%s -m > ./ci/output/rt_%s_%d.dat 2>&1", test_name, test_name, counter);
  system(command_buffer3);

  char command_buffer4[256];
  snprintf(command_buffer4, 256, "cd ../; mv ./%s-stat.json ci/output/%s-stat_%d.json", test_name, test_name, counter);
  system(command_buffer4);

  printf("Finished %s.\n\n", test_name_human);
}

void
analyzeTestOutput(const char* test_name)
{
  int counter = 0;

  char counter_buffer[64];
  snprintf(counter_buffer, 64, "output/%s_counter.dat", test_name);
  FILE *counter_ptr = fopen(counter_buffer, "r");
  if (counter_ptr != NULL) {
    fscanf(counter_ptr, "%d", &counter);
  }
  fclose(counter_ptr);

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

  for (int i = 1; i < counter + 1; i++) {
    char *output;
    long file_size;
    char buffer[128];
    snprintf(buffer, 128, "output/rt_%s_%d.dat", test_name, i);

    FILE *output_ptr = fopen(buffer, "rb");
    fseek(output_ptr, 0, SEEK_END);
    file_size = ftell(output_ptr);
    rewind(output_ptr);
    output = malloc(file_size * (sizeof(char)));
    fread(output, sizeof(char), file_size, output_ptr);
    fclose(output_ptr);

    updatecalls[i] = 0;
    if (strstr(output, "Number of update calls ") != NULL) {
      char *full_substring = strstr(output, "Number of update calls ");
      char substring[64];
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Number of update calls ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Number of update calls ")];
        substring_index += 1;
      }

      char *end_ptr;
      updatecalls[i] = strtol(substring, &end_ptr, 10);
    }

    forwardeuler[i] = 0;
    if (strstr(output, "Number of forward-Euler calls ") != NULL) {
      char *full_substring = strstr(output, "Number of forward-Euler calls ");
      char substring[64];
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Number of forward-Euler calls ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Number of forward-Euler calls ")];
        substring_index += 1;
      }

      char *end_ptr;
      forwardeuler[i] = strtol(substring, &end_ptr, 10);
    }

    rk2failures[i] = 0;
    if (strstr(output, "Number of RK stage-2 failures ") != NULL) {
      char *full_substring = strstr(output, "Number of RK stage-2 failures ");
      char substring[64];
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Number of RK stage-2 failures ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Number of RK stage-2 failures ")];
        substring_index += 1;
      }

      char *end_ptr;
      rk2failures[i] = strtol(substring, &end_ptr, 10);
    }

    rk3failures[i] = 0;
    if (strstr(output, "Number of RK stage-3 failures ") != NULL) {
      char *full_substring = strstr(output, "Number of RK stage-3 failures ");
      char substring[64];
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Number of RK stage-3 failures ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Number of RK stage-3 failures ")];
        substring_index += 1;
      }

      char *end_ptr;
      rk3failures[i] = strtol(substring, &end_ptr, 10);
    }

    speciesrhs[i] = 0.0;
    if (strstr(output, "Species RHS calc took ") != NULL) {
      char *full_substring = strstr(output, "Species RHS calc took ");
      char substring[64];
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Species RHS calc took ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Species RHS calc took ")];
        substring_index += 1;
      }

      char *end_ptr;
      speciesrhs[i] = strtod(substring, &end_ptr);
    }

    speciescollisionsrhs[i] = 0.0;
    if (strstr(output, "Species collisions RHS calc took ") != NULL) {
      char *full_substring = strstr(output, "Species collisions RHS calc took ");
      char substring[64];
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Species collisions RHS calc took ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Species collisions RHS calc took ")];
        substring_index += 1;
      }

      char *end_ptr;
      speciescollisionsrhs[i] = strtod(substring, &end_ptr);
    }

    fieldrhs[i] = 0.0;
    if (strstr(output, "Field RHS calc took ") != NULL) {
      char *full_substring = strstr(output, "Field RHS calc took ");
      char substring[64];
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Field RHS calc took ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Field RHS calc took ")];
        substring_index += 1;
      }

      char *end_ptr;
      fieldrhs[i] = strtod(substring, &end_ptr);
    }

    speciescollisionalmoments[i] = 0.0;
    if (strstr(output, "Species collisional moments took ") != NULL) {
      char *full_substring = strstr(output, "Species collisional moments took ");
      char substring[64];
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Species collisional moments took ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Species collisional moments took ")];
        substring_index += 1;
      }

      char *end_ptr;
      speciescollisionalmoments[i] = strtod(substring, &end_ptr);
    }

    totalupdate[i] = 0.0;
    if (strstr(output, "Total updates took ") != NULL) {
      char *full_substring = strstr(output, "Total updates took ");
      char substring[64];
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Total updates took ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Total updates took ")];
        substring_index += 1;
      }

      char *end_ptr;
      totalupdate[i] = strtod(substring, &end_ptr);
    }
    
    char *temp = output;
    memoryleakcount[i] = 0;
    memoryleaks[i] = (char*)malloc(1024 * sizeof(char));
    while (strstr(temp, "0x") != NULL) {
      temp = strstr(temp, "0x");

      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }

      int substring_index = 0;
      while (temp[substring_index] != ' ' && temp[substring_index] != '\n') {
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
      if (count == 1) {
        memoryleakcount[i] += 1;
        memoryleaks[i] = strcat(memoryleaks[i], substring);
        memoryleaks[i] = strcat(memoryleaks[i], " ");
      }
      
      temp += 1;
    }
  }

  for (int i = 1; i < counter + 1; i++) {
    printf("Build number: %d\n", i);
    if (i == 1) {
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
        printf("Species RHS time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", speciesrhs[i]);
      }
      else {
        printf("Species RHS time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", speciesrhs[i]);
      }

      if (speciescollisionsrhs[i] > speciescollisionsrhs[i - 1]) {
        printf("Species collision RHS time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", speciescollisionsrhs[i]);
      }
      else {
        printf("Species collision RHS time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", speciescollisionsrhs[i]);
      }

      if (fieldrhs[i] > fieldrhs[i - 1]) {
        printf("Field RHS time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", fieldrhs[i]);
      }
      else {
        printf("Field RHS time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", fieldrhs[i]);
      }

      if (speciescollisionalmoments[i] > speciescollisionalmoments[i - 1]) {
        printf("Species collisional moments time: " ANSI_COLOR_RED "%f" ANSI_COLOR_RESET "\n", speciescollisionalmoments[i]);
      }
      else {
        printf("Species collisional moments time: " ANSI_COLOR_GREEN "%f" ANSI_COLOR_RESET "\n", speciescollisionalmoments[i]);
      }

      if (memoryleakcount[i] != 0) {
        printf("Memory leaks: " ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", memoryleaks[i]);
      }
      else {
        printf("Memory leaks: " ANSI_COLOR_GREEN "None" ANSI_COLOR_RESET "\n");
      }

      if ((updatecalls[i] != updatecalls[i - 1]) || (forwardeuler[i] != forwardeuler[i - 1]) || (rk2failures[i] != rk2failures[i - 1]) || (rk3failures[i] != rk3failures[i - 1])) {
        printf("Correct: " ANSI_COLOR_RED "No" ANSI_COLOR_RESET "\n\n");
      }
      else {
        printf("Correct: " ANSI_COLOR_GREEN "Yes" ANSI_COLOR_RESET "\n\n");
      }
    }
  }
}

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
    runTest("gk_sheath_1x2v_p1", "1x2v Sheath Boundary Test with p = 1");
    runTest("gk_sheath_2x2v_p1", "2x2v Sheath Boundary Test with p = 1");
    runTest("gk_sheath_3x2v_p1", "3x2v Sheath Boundary Test with p = 1");
    runTest("gk_lapd_cart_3x2v_p1", "3x2v LAPD Test (in Cartesian coordinates) with p = 1");
    runTest("gk_lapd_cyl_3x2v_p1", "3x2v LAPD Test (in cylindrical coordinates) with p = 1");
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
      analyzeTestOutput("gk_sheath_1x2v_p1");
    }
    else if (option2 == 2) {
      printf("2x2v Sheath Boundary Test with p = 1:\n\n");
      analyzeTestOutput("gk_sheath_2x2v_p1");
    }
    else if (option2 == 3) {
      printf("3x2v Sheath Boundary Test with p = 1:\n\n");
      analyzeTestOutput("gk_sheath_3x2v_p1");
    }
    else if (option2 == 4) {
      printf("3x2v LAPD Test (in Cartesian coordinates) with p = 1:\n\n");
      analyzeTestOutput("gk_lapd_cart_3x2v_p1");
    }
    else if (option2 == 5) {
      printf("3x2v LAPD Test (in cylindrical coordinates) with p = 1:\n\n");
      analyzeTestOutput("gk_lapd_cyl_3x2v_p1");
    }
  }
  else if (option == 3) {
    system("rm -rf output");
    system("mkdir output");

    runTest("gk_sheath_1x2v_p1", "1x2v Sheath Boundary Test with p = 1");
    runTest("gk_sheath_2x2v_p1", "2x2v Sheath Boundary Test with p = 1");
    runTest("gk_sheath_3x2v_p1", "3x2v Sheath Boundary Test with p = 1");
    runTest("gk_lapd_cart_3x2v_p1", "3x2v LAPD Test (in Cartesian coordinates) with p = 1");
    runTest("gk_lapd_cyl_3x2v_p1", "3x2v LAPD Test (in cylindrical coordinates) with p = 1");
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
      runTest("gk_sheath_1x2v_p1", "1x2v Sheath Boundary Test with p = 1");
    }
    else if (option2 == 2) {
      runTest("gk_sheath_2x2v_p1", "2x2v Sheath Boundary Test with p = 1");
    }
    else if (option2 == 3) {
      runTest("gk_sheath_3x2v_p1", "3x2v Sheath Boundary Test with p = 1");
    }
    else if (option2 == 4) {
      runTest("gk_lapd_cart_3x2v_p1", "3x2v LAPD Test (in Cartesian coordinates) with p = 1");
    }
    else if (option2 == 5) {
      runTest("gk_lapd_cyl_3x2v_p1", "3x2v LAPD Test (in cylindrical coordinates) with p = 1");
    }
  }

  return 0;
}