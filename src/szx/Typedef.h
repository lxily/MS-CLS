////////////////////////////////
/// ���:	1.	�������ͱ���.
/// 
/// ��ע:	1.	
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

// ��ʶ��. �� 0 ��ʼ������ŵ�����.
using ID = int;

// ���� (�������ʱ�䵥λ). ������ȵ���Сʱ�䵥λ.
using Step = int;

// ����. �źŴ������С��������λ.
using Bit = long long;

// ���� (��ʵ����ʱ�䵥λ). �㷨����ʱ��. ���ڳ�����ֹ�жϻ���־��¼��.
using Millisecond = long long;

// ������������������ֵ.
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
