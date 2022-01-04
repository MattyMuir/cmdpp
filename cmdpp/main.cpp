#include <mpir.h>
#include <primesieve.hpp>
#include <iostream>
#include <thread>
#include <sstream>

#include "mpzArray.h"
#include "Timer.h"

#define CONSTANT_WIDTH_BATCHES 0
#define LOG_STATUS 1
#define FALLBACK_VALUES 1

struct Range
{
    Range() {}
    Range(uint32_t start_, uint32_t end_)
        : start(start_), end(end_) {}
    uint32_t start = 0, end = 0;
};

void BatchRanges(Range splitRange, int n, std::vector<Range>& ranges)
{
    ranges = std::vector<Range>(n);
    long double w = splitRange.end - splitRange.start + 1;

#if CONSTANT_WIDTH_BATCHES == 1
    int batchSize = w / n;
    for (int t = 0; t < n; t++)
    {
        ranges[t].start = splitRange.start + batchSize * t;
        ranges[t].end = splitRange.start + ((t < n - 1) ? (batchSize * (t + 1) - 1) : w - 1);
    }
#else
    long double a = 1;

    for (int i = 0; i < n - 1; i++)
    {
        a = (a + sqrt(a * a + 4)) / 2;
    }

    long double x1 = w / a;

    a = 0;
    for (int i = 0; i < n; i++)
    {

        long double start = x1 * a;
        a = (a + sqrt(a * a + 4)) / 2;
        long double end = x1 * a;

        ranges[i].start = splitRange.start + ceil(start);
        ranges[i].end = splitRange.start + ((i == n - 1) ? (w - 1) : (floor(end)));
    }
#endif
}

void CheckBatch(Range check, int arrOffset, mpz_t& startProd, std::vector<uint64_t>& primes)
{
#if LOG_STATUS == 1
    Timer t;
#endif
    for (int index = check.start; index <= check.end; index++)
    {
        mpz_mul_ui(startProd, startProd, primes[index - 1 - arrOffset]);
        mpz_add_ui(startProd, startProd, 1);
        if (mpz_divisible_ui_p(startProd, primes[index - arrOffset]))
        {
            std::string nStr;
            nStr = std::to_string(index) + "\n";
            std::cout << nStr;
        }
        mpz_sub_ui(startProd, startProd, 1);
    }

#if LOG_STATUS == 1
    t.Stop();
    std::stringstream msg;
    msg << "Thread: " << check.start << " - " << check.end << " finished in " << t.duration * 0.000001 << "s" << std::endl;
    std::cout << msg.str();
#endif
}

void Run(Range r, std::vector<Range>* rangesPtr = nullptr)
{
    uint32_t N = r.end - r.start + 1;

    // Hardware lookup
    int threadNum = std::thread::hardware_concurrency();

    Timer t;

    // Generate end + 1 primes
    std::vector<uint64_t> prePrimes;
    std::vector<uint64_t> primes;

    prePrimes.reserve(r.start - 1);
    primes.reserve(N + 1);
    primesieve::iterator it;
    for (int p = 0; p < r.start - 1; p++)
    {
        prePrimes.push_back(it.next_prime());
    }
    for (int p = r.start - 1; p < r.end + 1; p++)
    {
        primes.push_back(it.next_prime());
    }

    // Copy N primes into arr
    mpzArray preArr;
    mpzArray arr;

    preArr.Reserve(r.start - 1);
    arr.Reserve(N);
    for (uint64_t p : prePrimes) { preArr.BackUI(p); }
    for (int p = 0; p < primes.size() - 1; p++) { arr.BackUI(primes[p]); }
    prePrimes.clear();

    // Calculate start offset product
    int activeNum = preArr.Size();

    while (activeNum > 1)
    {
        for (int i = 0; i < activeNum / 2; i++)
        {
            mpz_mul(preArr[i], preArr[i], preArr[activeNum - i - 1]);
        }

        activeNum = (activeNum + 1) / 2;
        preArr.Resize(activeNum);
    }
    mpz_t offsetProd;
    if (preArr.Size() > 0)
    {
        mpz_init_set(offsetProd, preArr[0]);
    }
    else
    {
        mpz_init_set_ui(offsetProd, 1);
    }
    preArr.Clear();
    std::cout << "Offset Product Calculated" << std::endl;

    // Pre-calculate start products
    // Calculate batch ranges
    std::vector<Range> ranges;
    if (rangesPtr)
    {
        ranges = *rangesPtr;
    }
    else
    {
        BatchRanges(r, threadNum, ranges);
    }

    // Loop through batches
    for (int b = 0; b < threadNum - 1; b++)
    {
        int startIndex = ranges[b].start - r.start;
        int endIndex = ranges[b].end - r.start;

        // Multiply batch in arr
        activeNum = endIndex - startIndex + 1;

        while (activeNum > 1)
        {
            for (int i = 0; i < activeNum / 2; i++)
            {
                mpz_mul(arr[startIndex + i], arr[startIndex + i], arr[startIndex + activeNum - i - 1]);
            }

            activeNum = (activeNum + 1) / 2;
        }

        if (b != 0)
        {
            mpz_mul(arr[startIndex], arr[startIndex], arr[ranges[b - 1].start - r.start]);
        }
    }

    // Multiply products by offset
    for (int b = 0; b < threadNum - 1; b++)
    {
        mpz_mul(arr[ranges[b].start - r.start], arr[ranges[b].start - r.start], offsetProd);

        // Shift final start products towards the first 'threadNum - 1' spaces in 'arr'
        if (b > 0)
        {
            mpz_init_set(arr[b], arr[ranges[b].start - r.start]);
        }
    }

    // Trim off rest of the array only leaving the calculated start products
    arr.Resize(threadNum - 1);

    // Offset batches by offset prod
    std::cout << "Start Points Calculated" << std::endl;

    // Initialize threads on each batch
    std::vector<std::thread> threads;

    for (int t = 0; t < threadNum; t++)
    {
        int startIndex = ranges[t].start;
        int endIndex = ranges[t].end;

        if (t == 0)
        {
            threads.push_back(std::thread(CheckBatch, ranges[t], r.start - 1, std::ref(offsetProd), std::ref(primes)));
        }
        else
        {
            threads.push_back(std::thread(CheckBatch, ranges[t], r.start - 1, std::ref(arr[t - 1]), std::ref(primes)));
        }
    }

    // Wait till all threads are complete
    for (std::thread& t : threads) { t.join(); }

    t.Stop();
    std::cout << "Took: " << t.duration * 0.001 << "ms" << std::endl;
}

/// Generate a way of splitting 'r' into 'n' ranges with similar computation times
std::vector<Range> Tune(Range r, int n, int iter)
{

}

int main(int argc, char* argv[])
{
    uint32_t start, end;

    if (argc != 3)
    {
#if FALLBACK_VALUES == 0
        std::cout << "Wrong number of args.\n";
#else
        start = 1; end = 1e5;
        goto run;
#endif
        return 1;
    }
    
    try
    {
        start = std::stoi(argv[1]);
        end = std::stoi(argv[2]);
    }
    catch (const std::invalid_argument& e)
    {
        std::cout << "Invalid argument.\n";
        return 1;
    }
    catch (const std::out_of_range& e)
    {
        std::cout << "Argument out of range.\n";
        return 1;
    }

run:
    Tune(Range(start, end), std::thread::hardware_concurrency(), 10);

    //Run(Range(start, end));
}