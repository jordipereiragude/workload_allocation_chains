/*
 * Copyright (c) 2023 Jordi Pereira, Marcus Ritt
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#include "memory.hpp"

#include <fstream>
using namespace std;

#ifdef __linux__
#include <sys/types.h>
#include <sys/sysinfo.h>
#include <unistd.h>
#include <linux/sysctl.h>
#endif

#ifdef __APPLE__
#include <mach/task.h>
#include <mach/mach_init.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif

#ifdef _WINDOWS
#include <windows.h>
#else
#include <sys/resource.h>
#endif

size_t usedVMproc(bool resident) {
#if defined(__linux__)
  size_t size = 0;
  ifstream file("/proc/self/statm");
  if (file.good()) {
    file >> size;
    size *= getpagesize();
  }
  return size;
#elif defined(__APPLE__)
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
  task_info(current_task(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count);
  size_t size = (resident ? t_info.resident_size : t_info.virtual_size);
  return size;
#elif defined(_WINDOWS)
  PROCESS_MEMORY_COUNTERS counters;
  if (GetProcessMemoryInfo (GetCurrentProcess(), &counters, sizeof (counters)))
    return counters.PagefileUsage;
  else return 0;
#else
  return 0;
#endif
}

long long totalVM() {
#if defined(__linux__)
  struct sysinfo memInfo;
  sysinfo (&memInfo);

  long long totalVirtualMem = memInfo.totalram;
  //Add other values in next statement to avoid int overflow on right hand side...
  totalVirtualMem += memInfo.totalswap;
  totalVirtualMem *= memInfo.mem_unit;
  return totalVirtualMem;
#else
  return 0;
#endif
}

long long usedVM() {
#if defined(__linux__)
  struct sysinfo memInfo;
  sysinfo (&memInfo);

  long long virtualMemUsed = memInfo.totalram - memInfo.freeram;
  virtualMemUsed += memInfo.totalswap - memInfo.freeswap;
  virtualMemUsed *= memInfo.mem_unit;
  return virtualMemUsed;
#elif defined(__APPLE__)
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  return((long long)(usage.ru_maxrss));
#else
  return 0;
#endif

}

int64_t totalM() {
#if defined(DOESNT_WORK)
  int mib[2];
  int64_t physical_memory;
  mib[0] = CTL_HW;
  mib[1] = HW_MEMSIZE;
  length = sizeof(int64_t);
  sysctl(mib, 2, &physical_memory, &length, NULL, 0);
  return physical_memory;
#else
  ifstream file("/proc/meminfo");
  string dummy;
  int64_t mem;
  file >> dummy >> mem;
  return mem/1024;
#endif
}
