////////////////////////////////
/// ¼ò½é:	1.	¼ÆÊ±Æ÷.
/// 
/// ±¸×¢:	1.	
////////////////////////////////

#ifndef CN_HUST_HAFV_UTIL_TIMER_H
#define CN_HUST_HAFV_UTIL_TIMER_H


#include <chrono>
#include <ctime>


namespace hust {
namespace util {

class TimerBase {
public:
    using Microsecond = unsigned long long;


    static constexpr Microsecond MicrosecondsPerSecond = 1000000;
    static constexpr double MillisecondsPerSecond = 1000;
    static constexpr double ClocksPerSecond = CLOCKS_PER_SEC;
    static constexpr clock_t ClocksPerMillisecond = CLOCKS_PER_SEC / static_cast<clock_t>(MillisecondsPerSecond);


    // there is no need to free the pointer. the format of the format string is 
    // the same as std::strftime() in http://en.cppreference.com/w/cpp/chrono/c/strftime.
    static const char* getLocalTime(const char *format = "%Y-%m-%d(%a)%H:%M:%S") {
        static constexpr int DateBufSize = 64;
        static char buf[DateBufSize];
        time_t t = time(NULL);
        tm *date = localtime(&t);
        strftime(buf, DateBufSize, format, date);
        return buf;
    }
    static const char* getTightLocalTime() { return getLocalTime("%Y%m%d%H%M%S"); }
};

class TimerCpp : public TimerBase {
public:
    using Millisecond = std::chrono::milliseconds;
    using TimePoint = std::chrono::steady_clock::time_point;
    using Clock = std::chrono::steady_clock;


    TimerCpp(const Millisecond &duration = Millisecond(0), const TimePoint &st = Clock::now())
        : startTime(st), endTime(startTime + duration) {
    }


    static Millisecond toMillisecond(double second) {
        return Millisecond(static_cast<long long>(second * MillisecondsPerSecond));
    }

    static Millisecond durationInMillisecond(const TimePoint &start, const TimePoint &end) {
        return std::chrono::duration_cast<Millisecond>(end - start);
    }

    static double durationInSecond(const TimePoint &start, const TimePoint &end) {
        return static_cast<double>(std::chrono::duration_cast<Millisecond>(end - start).count()) / MillisecondsPerSecond;
    }

    bool isTimeOut() const {
        return (Clock::now() > endTime);
    }

    Millisecond restMilliseconds() const {
        return durationInMillisecond(Clock::now(), endTime);
    }

    double restSeconds() const {
        return durationInSecond(Clock::now(), endTime);
    }

    Millisecond elapsedMilliseconds() const {
        return durationInMillisecond(startTime, Clock::now());
    }

    double elapsedSeconds() const {
        return durationInSecond(startTime, Clock::now());
    }

    const TimePoint& getStartTime() const { return startTime; }
    const TimePoint& getEndTime() const { return endTime; }

protected:
    TimePoint startTime;
    TimePoint endTime;
};

class TimerC : public TimerBase {
public:
    using Millisecond = clock_t;
    using TimePoint = clock_t;
    struct Clock { static TimePoint now() { return clock(); } };


    TimerC(const Millisecond &duration = Millisecond(0), const TimePoint &st = Clock::now())
        : startTime(st), endTime(startTime + duration * ClocksPerMillisecond) {
    }


    static Millisecond toMillisecond(double second) {
        return static_cast<Millisecond>(second * MillisecondsPerSecond);
    }

    static Millisecond durationInMillisecond(const TimePoint &start, const TimePoint &end) {
        return (end - start) / ClocksPerMillisecond;
    }

    static double durationInSecond(const TimePoint &start, const TimePoint &end) {
        return static_cast<double>(end - start) / ClocksPerSecond;
    }

    bool isTimeOut() const {
        return (Clock::now() > endTime);
    }

    Millisecond restMilliseconds() const {
        return durationInMillisecond(Clock::now(), endTime);
    }

    double restSeconds() const {
        return durationInSecond(Clock::now(), endTime);
    }

    Millisecond elapsedMilliseconds() const {
        return durationInMillisecond(startTime, Clock::now());
    }

    double elapsedSeconds() const {
        return durationInSecond(startTime, Clock::now());
    }

    const TimePoint& getStartTime() const { return startTime; }
    const TimePoint& getEndTime() const { return endTime; }

protected:
    TimePoint startTime;
    TimePoint endTime;
};


using Timer = TimerCpp;

}
}


#endif // CN_HUST_HAFV_UTIL_TIMER_H
