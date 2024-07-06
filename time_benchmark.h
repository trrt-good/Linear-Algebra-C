#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#if defined(_WIN32) || defined(_WIN64)
    #include <windows.h>
#elif defined(__APPLE__)

    #include <mach/mach_time.h>
#else 
    #error "OS not supported"
#endif

struct time_tracker{
    #if defined(_WIN32) || defined(_WIN64)
        LARGE_INTEGER tps;
        LARGE_INTEGER t1, t2;
    #elif defined(__APPLE__)
        mach_timebase_info_data_t timebaseInfo;
        uint64_t t1, t2;
    #endif
    float timeDiff;
};

struct time_tracker* tracker;

// init the time tracker struct based on the OS
void init_time_tracker(){
    srand(time(NULL));
    struct time_tracker* tt = (struct time_tracker*)malloc(sizeof(struct time_tracker));
    #if defined(_WIN32) || defined(_WIN64)
        QueryPerformanceFrequency(&tt->tps);
    #elif defined(__APPLE__)
        mach_timebase_info(&tt->timebaseInfo);
    #endif 
    tracker = tt;
};

// start tracking time
void start_tracking_time(){
    init_time_tracker();
    #if defined(_WIN32) || defined(_WIN64)
        QueryPerformanceCounter(&tracker->t1);
    #elif defined(__APPLE__)
        tracker->t1 = mach_absolute_time();
    #endif
}

// end tracking time
void end_tracking_time(){
    #if defined(_WIN32) || defined(_WIN64)
        QueryPerformanceCounter(&tracker->t2);
        tracker->timeDiff = (tracker->t2.QuadPart - tracker->t1.QuadPart) * 1000.0 / tracker->tps.QuadPart;
    #elif defined(__APPLE__)
        tracker->t2 = mach_absolute_time();
        uint64_t elapsed = tracker->t2 - tracker->t1;
        tracker->timeDiff = (double)elapsed * tracker->timebaseInfo.numer / tracker->timebaseInfo.denom / 1e6; // Convert to milliseconds
    #endif
}

// get the time taken in milliseconds
float get_time_taken(){
    return tracker->timeDiff;
}