#pragma once

// Status of update() method. If success is 0 (failure) then the
// simulation needs to be aborted and can't continue.
struct gkyl_update_status {
    int success; // 1 if update worked, 0 if a fatal error
    double dt_actual; // actual time-step taken
    double dt_suggested; // suggested stable time-step
};

