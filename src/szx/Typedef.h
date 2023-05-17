////////////////////////////////
/// 简介:	1.	常用类型别名.
/// 
/// 备注:	1.	
////////////////////////////////

#ifndef CN_HUST_HAFV_UTIL_TYPEDEF_H
#define CN_HUST_HAFV_UTIL_TYPEDEF_H


#include <string>
#include <array>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <functional>
#include <cstdint>

#include "Arr2D.h"


namespace hust {

// 标识符. 从 0 开始连续编号的整数.
using ID = int;

// 节拍 (仿真过程时间单位). 任务调度的最小时间单位.
using Step = int;

// 比特. 信号传输的最小数据量单位.
using Bit = long long;

// 毫秒 (现实世界时间单位). 算法运行时间. 用于程序终止判断或日志记录等.
using Millisecond = long long;

// 计数或其他无量纲数值.
using Int = long long;
using Real = double;

using Str = std::string;
using StrPos = Str::size_type;

template<typename T, size_t Size>
using Arr = std::array<T, Size>;

template<typename T, typename IndexType = int>
using Arr2D = util::Array2D<T, IndexType>;

template<typename T>
using Vec = std::vector<T>;

template<typename T>
using Set = std::set<T>;

template<typename Key, typename Value>
using Map = std::map<Key, Value>;

template<typename Key, typename Value>
using HashMap = std::unordered_map<Key, Value>;


template<typename Prototype>
using Func = std::function<Prototype>;


struct Void {};


template<typename T>
const char* toBytes(const T *ptr) { return reinterpret_cast<const char*>(ptr); }

template<typename T>
char* toBytes(T *ptr) { return reinterpret_cast<char*>(ptr); }


static constexpr ID InvalidId = -1;

}


#endif // CN_HUST_HAFV_UTIL_TYPEDEF_H
