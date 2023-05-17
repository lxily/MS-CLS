////////////////////////////////
/// usage : 1.	configurable priority queue.
/// 
/// note  : 1.	for all implementations marked with [ValueAsKey] instead of [SeperatedKey],
///				the users must not alter the behavior of the `LessKeyPred` between each pair of items.
/// 			e.g., if `lessPred = [](int l, int r) { return a[l] < a[r]; }` then `++a[i]` is forbidden
/// 			unless item `i` is removed from the queue.
///				an alternative is to explicitly store each key in its pairing value and give up updating the keys.
/// 		2.	for all implementations that do not support removing items or altering the keys,
/// 			there must be additional verification mechanisms to detect potentially out-of-date top-priority keys.
/// 		3.	Additional check for calling `top()` or `pop()` on empty queue is not implemented.
/// 			Users should call `empty()` manually or runtime errors will occur.
////////////////////////////////

#ifndef CN_HUST_HAFV_UTIL_PRIORITY_QUEUE_H
#define CN_HUST_HAFV_UTIL_PRIORITY_QUEUE_H


#include <algorithm>
#include <functional>
#include <vector>
#include <set>
#include <map>
#include <type_traits>

#include <HAFV/Util/Flag.h>
#include <HAFV/Util/Typedef.h>
#include <HAFV/Util/Math.h>


namespace hust {
namespace hafv {

// smaller key means higher priority.
// it is effective on radix based implementations only and will not bother compare based ones.
using Priority = int;

template<typename T>
using LessKeyCmp = std::function<bool(const T &l, const T &r)>;

template<typename T>
using GreaterKeyCmp = std::function<bool(const T &l, const T &r)>;


// [CompareBased][ValueAsKey][NoDup][IncreaseKey][DecreaseKey]
template<typename T, typename LessKeyPred = LessKeyCmp<T>>
class PriorityQueueBySet {
public:
    using ContainerType = std::set<T, LessKeyPred>;


    static constexpr Priority DefaultHintItemNum = 64;


    explicit PriorityQueueBySet(const LessKeyPred &lessPred)
        : container(lessPred) {
    }

    explicit PriorityQueueBySet(Priority /*hintItemNum*/ = DefaultHintItemNum, const LessKeyPred &lessPred = std::less<T>())
        : container(lessPred) {
    }

    explicit PriorityQueueBySet(const ContainerType &itemContainer, const LessKeyPred &lessPred = std::less<T>())
        : container(itemContainer.begin(), itemContainer.end(), lessPred) {
    }


// peek the item with the minimal key (highest priority).
    T top() const { return *container.begin(); }

    Priority topPriority() const {
        static_assert(std::is_convertible<T, Priority>::value, "Please implement a converter from the value type to the key type.");
        return static_cast<Priority>(top());
    }

    // remove the item with the minimal key (highest priority).
    void pop() { container.erase(container.begin()); }

    // put the item in proper place regarding its precedence.
    void push(const T &e, Priority /*key*/) { container.insert(e); }

    template<typename RandomNumGenerator>
    void push(const T &e, Priority key, RandomNumGenerator & /*rgen*/) { push(e, key); }

    // alter the key of an existing item.
    void update(const T &e, Priority key) {
        remove(e);
        push(e, key);
    }
    void update(const T &e, Priority key, std::function<void(void)> alterLessKeyPredBehavior) {
        remove(e);
        alterLessKeyPredBehavior();
        push(e, key);
    }

    // make the key of an existing item smaller (higher priority).
    void descend(const T &e, Priority key) { update(e, key); }

    // make the key of an existing item larger (lower priority).
    void ascend(const T &e, Priority key) { update(e, key); }

    // remove the given item.
    void remove(const T &e) { container.erase(e); }

    // reorganize to meet the requirement of the priority queue.
    void reorder() {} // nothing to be done.

    // pop all items.
    void clear() { container.clear(); }

    // if there is no item in the queue.
    bool empty() const { return container.empty(); }

    Priority size() const { return static_cast<Priority>(container.size()); }

    void reserve(Priority /*capacity*/, Priority /*key*/) {}
    void reserve(Priority /*capacity*/) {}

    const ContainerType& getContainer() const { return container; }

protected:
    ContainerType container;
};


// [CompareBased][SeperatedKey]
template<typename T, typename LessKeyPred = LessKeyCmp<T>>
class PriorityQueueByMap {
public:
    using ContainerType = std::map<Priority, Vec<T>, LessKeyPred>;
    // EXTEND[szx][8]: use `std::multimap` or replace `Vec` with `std::set` to enable remove operation?


    static constexpr Priority DefaultHintItemNum = 64;


    explicit PriorityQueueByMap(const LessKeyPred &lessPred)
        : container(lessPred), count(0) {
    }

    explicit PriorityQueueByMap(Priority /*hintItemNum*/ = DefaultHintItemNum, const LessKeyPred &lessPred = std::less<T>())
        : container(lessPred), count(0) {
    }

    explicit PriorityQueueByMap(const ContainerType &itemContainer, const LessKeyPred &lessPred = std::less<T>())
        : container(itemContainer.begin(), itemContainer.end(), lessPred) {
    }


// peek the item with the minimal key (highest priority).
    T top() const { return container.begin()->second.back(); }

    Priority topPriority() const { return container.begin()->first; }

    // remove the item with the minimal key (highest priority).
    void pop() {
        container.begin()->second.pop_back();
        if (container.begin()->second.empty()) { container.erase(container.begin()); }
        --count;
    }

    // put the item in proper place regarding its precedence.
    void push(const T &e, Priority key) {
        container[key].push_back(e);
        ++count;
    }

    template<typename RandomNumGenerator>
    void push(const T &e, Priority key, RandomNumGenerator &rgen) {
        if (container[key].empty()) { push(e, key); return; }
        size_t i = static_cast<size_t>(rgen()) % container[key].size();
        container[key].push_back(container[key][i]);
        container[key][i] = e;
        ++count;
    }

    // alter the key of an existing item.
    void update(const T & /*e*/, Priority /*key*/); // { throw NotImplementedException(); }
    void update(const T & /*e*/, Priority /*key*/, std::function<void(void)> /*alterLessKeyPredBehavior*/); // { throw NotImplementedException(); }

    // make the key of an existing item smaller (higher priority).
    void descend(const T & /*e*/, Priority /*key*/); // { throw NotImplementedException(); }

    // make the key of an existing item larger (lower priority).
    void ascend(const T & /*e*/, Priority /*key*/); // { throw NotImplementedException(); }

    // remove the given item.
    void remove(const T & /*e*/); // { throw NotImplementedException(); }

    // reorganize to meet the requirement of the priority queue.
    void reorder() {} // nothing to be done.

    // pop all items.
    void clear() {
        container.clear();
        count = 0;
    }

    // if there is no item in the queue.
    bool empty() const { return count == 0; }

    Priority size() const { return count; }

    void reserve(Priority capacity, Priority key) { container[key].reserve(capacity); }
    void reserve(Priority /*capacity*/) {}

    const ContainerType& getContainer() const { return container; }

protected:
    ContainerType container;
    // the number of items.
    Priority count;
};


// EXT[szx][4]: the template parameter and the constructor parameter are different from other implementations.
//              the consistency needs to be improved.
// [CompareBased][ValueAsKey]
template<typename T, typename GreaterKeyPred = GreaterKeyCmp<T>>
class PriorityQueueByHeap {
public:
    using ContainerType = Vec<T>;


    static constexpr Priority DefaultHintItemNum = 64;


    // use `greater()` as the comparer to achieve a min heap (STL only provide max heap).
    explicit PriorityQueueByHeap(const GreaterKeyPred &greaterPred)
        : greater(greaterPred) {
    }

    explicit PriorityQueueByHeap(Priority hintItemNum = DefaultHintItemNum, const GreaterKeyPred &greaterPred = std::greater<T>())
        : greater(greaterPred) {
        container.reserve(hintItemNum);
    }

    explicit PriorityQueueByHeap(const ContainerType &itemContainer, const GreaterKeyPred &greaterPred = std::greater<T>())
        : container(itemContainer), greater(greaterPred) {
    }


// peek the item with the minimal key (highest priority).
    T top() const { return container.front(); }

    Priority topPriority() const {
        static_assert(std::is_convertible<T, Priority>::value, "Please implement a converter from the value type to the key type.");
        return static_cast<Priority>(top());
    }

    // remove the item with the minimal key (highest priority).
    void pop() {
        std::pop_heap(container.begin(), container.end(), greater);
        container.pop_back();
    }

    // put the item in proper place regarding its precedence.
    void push(const T &e, Priority /*key*/) {
        container.push_back(e);
        std::push_heap(container.begin(), container.end(), greater);
    }

    template<typename RandomNumGenerator>
    void push(const T &e, Priority key, RandomNumGenerator & /*rgen*/) { push(e, key); }

    // alter the key of an existing item.
    void update(const T & /*e*/, Priority /*key*/); // { throw NotImplementedException(); }
    void update(const T & /*e*/, Priority /*key*/, std::function<void(void)> /*alterLessKeyPredBehavior*/); // { throw NotImplementedException(); }

    // make the key of an existing item smaller (higher priority).
    void descend(const T & /*e*/, Priority /*key*/); // { throw NotImplementedException(); }

    // make the key of an existing item larger (lower priority).
    void ascend(const T & /*e*/, Priority /*key*/); // { throw NotImplementedException(); }

    // remove the given item.
    void remove(const T & /*e*/); // { throw NotImplementedException(); } // EXT[szx][8]: find and incrementally reorder?

    // reorganize to meet the requirement of the priority queue.
    void reorder() { std::make_heap(container.begin(), container.end(), greater); }

    // pop all items.
    void clear() { container.clear(); }

    // if there is no item in the queue.
    bool empty() const { return container.empty(); }

    Priority size() const { return static_cast<Priority>(container.size()); }

    void reserve(Priority /*capacity*/, Priority /*key*/) {}
    void reserve(Priority capacity) { container.reserve(capacity); }

    const ContainerType& getContainer() const { return container; }

protected:
    ContainerType container;
    GreaterKeyPred greater;
};


// [RadixBased][SeperatedKey][PositivePriorityOnly]
template<typename T>
class PriorityQueueByBucketL1FixedSizeMinTop {
public:
    using ContainerType = Vec<Vec<T>>;


    static constexpr Priority InvalidIndex = (std::numeric_limits<Priority>::max)() / 2;
    static constexpr Priority DefaultHintItemNum = 1024;


    explicit PriorityQueueByBucketL1FixedSizeMinTop(Priority maxBucketNum) : container(maxBucketNum),
        firstNonEmptyIndex(InvalidIndex), lastNonEmptyIndex(0), count(0) {
    }

    explicit PriorityQueueByBucketL1FixedSizeMinTop(const ContainerType &itemContainer) : container(itemContainer),
        firstNonEmptyIndex(InvalidIndex), lastNonEmptyIndex(0), count(0) {
    }


    // peek the item with the minimal key (highest priority).
    T top() { return container[updateFirstNonEmptyIndex()].back(); }

    Priority topPriority() { return updateFirstNonEmptyIndex(); }

    // remove the item with the minimal key (highest priority).
    void pop() {
        container[updateFirstNonEmptyIndex()].pop_back();
        --count;
    }

    // TODO[szx][9]: what if the priority is negative?
    //               pass an offset in constructor and add it to the priority in each push?
    // put the item in proper place regarding its precedence.
    void push(const T &e, Priority key) {
        container[key].push_back(e);
        if (key < firstNonEmptyIndex) { firstNonEmptyIndex = key; }
        if (key > lastNonEmptyIndex) { lastNonEmptyIndex = key; }
        ++count;
    }

    template<typename RandomNumGenerator>
    void push(const T &e, Priority key, RandomNumGenerator &rgen) {
        if (container[key].empty()) { push(e, key); return; }
        size_t i = static_cast<size_t>(rgen()) % container[key].size();
        container[key].push_back(container[key][i]);
        container[key][i] = e;
        if (key < firstNonEmptyIndex) { firstNonEmptyIndex = key; }
        if (key > lastNonEmptyIndex) { lastNonEmptyIndex = key; }
        ++count;
    }

    // alter the key of an existing item.
    void update(const T & /*e*/, Priority /*key*/); // { throw NotImplementedException(); }
    void update(const T & /*e*/, Priority /*key*/, std::function<void(void)> /*alterLessKeyPredBehavior*/); // { throw NotImplementedException(); }

    // make the key of an existing item smaller (higher priority).
    void descend(const T & /*e*/, Priority /*key*/); // { throw NotImplementedException(); }

    // make the key of an existing item larger (lower priority).
    void ascend(const T & /*e*/, Priority /*key*/); // { throw NotImplementedException(); }

    // remove the given item.
    void remove(const T & /*e*/); // { throw NotImplementedException(); }

    // reorganize to meet the requirement of the priority queue.
    void reorder() {} // EXT[szx][9]: clean invalid items.

    // pop all items.
    void clear() {
        for (Priority i = firstNonEmptyIndex; i <= lastNonEmptyIndex; ++i) { container[i].clear(); }
        firstNonEmptyIndex = InvalidIndex;
        lastNonEmptyIndex = 0;
        count = 0;
    }

    // if there is no item in the queue.
    bool empty() { return count == 0; }

    Priority size() const { return count; }

    const ContainerType& getContainer() const { return container; }
    const ContainerType& getContainer(Priority bucket) const { return container[bucket]; }

    void reserve(Priority bucket, int capacity) { container[bucket].reserve(capacity); }

protected:
    Priority updateFirstNonEmptyIndex() {
        for (; firstNonEmptyIndex <= lastNonEmptyIndex; ++firstNonEmptyIndex) {
            if (!container[firstNonEmptyIndex].empty()) { return firstNonEmptyIndex; }
        }
        return InvalidIndex;
    }


    // `container[bucket][i]` is the i_th item in the bucket.
    ContainerType container;
    // the index of the first non-empty group.
    Priority firstNonEmptyIndex;
    // the index of the last non-empty group.
    Priority lastNonEmptyIndex;
    // the number of items.
    Priority count;
};
// [RadixBased][SeperatedKey][PositivePriorityOnly]
template<typename T>
class PriorityQueueByBucketL1FixedSizeMaxTop {
public:
    using ContainerType = Vec<Vec<T>>;


    static constexpr Priority InvalidIndex = -1;
    static constexpr Priority DefaultHintItemNum = 1024;


    explicit PriorityQueueByBucketL1FixedSizeMaxTop(Priority maxBucketNum) : container(maxBucketNum),
        firstNonEmptyIndex(InvalidIndex), lastNonEmptyIndex(0), count(0) {
    }

    explicit PriorityQueueByBucketL1FixedSizeMaxTop(const ContainerType &itemContainer) : container(itemContainer),
        firstNonEmptyIndex(InvalidIndex), lastNonEmptyIndex(0), count(0) {
    }


    // peek the item with the minimal key (highest priority).
    T top() { return container[updateLastNonEmptyIndex()].back(); }

    Priority topPriority() { return updateLastNonEmptyIndex(); }

    // remove the item with the minimal key (highest priority).
    void pop() {
        container[updateLastNonEmptyIndex()].pop_back();
        --count;
    }

    // TODO[szx][9]: what if the priority is negative?
    //               pass an offset in constructor and add it to the priority in each push?
    // put the item in proper place regarding its precedence.
    void push(const T &e, Priority key) {
        container[key].push_back(e);
        if (key < firstNonEmptyIndex) { firstNonEmptyIndex = key; }
        if (key > lastNonEmptyIndex) { lastNonEmptyIndex = key; }
        ++count;
    }

    template<typename RandomNumGenerator>
    void push(const T &e, Priority key, RandomNumGenerator &rgen) {
        if (container[key].empty()) { push(e, key); return; }
        size_t i = static_cast<size_t>(rgen()) % container[key].size();
        container[key].push_back(container[key][i]);
        container[key][i] = e;
        if (key < firstNonEmptyIndex) { firstNonEmptyIndex = key; }
        if (key > lastNonEmptyIndex) { lastNonEmptyIndex = key; }
        ++count;
    }

    // alter the key of an existing item.
    void update(const T & /*e*/, Priority /*key*/); // { throw NotImplementedException(); }
    void update(const T & /*e*/, Priority /*key*/, std::function<void(void)> /*alterLessKeyPredBehavior*/); // { throw NotImplementedException(); }

    // make the key of an existing item smaller (higher priority).
    void descend(const T & /*e*/, Priority /*key*/); // { throw NotImplementedException(); }

    // make the key of an existing item larger (lower priority).
    void ascend(const T & /*e*/, Priority /*key*/); // { throw NotImplementedException(); }

    // remove the given item.
    void remove(const T & /*e*/); // { throw NotImplementedException(); }

    // reorganize to meet the requirement of the priority queue.
    void reorder() {} // EXT[szx][9]: clean invalid items.

    // pop all items.
    void clear() {
        for (Priority i = firstNonEmptyIndex; i <= lastNonEmptyIndex; ++i) { container[i].clear(); }
        firstNonEmptyIndex = InvalidIndex;
        lastNonEmptyIndex = 0;
        count = 0;
    }

    // if there is no item in the queue.
    bool empty() { return count == 0; }

    Priority size() const { return count; }

    const ContainerType& getContainer() const { return container; }
    const ContainerType& getContainer(Priority bucket) const { return container[bucket]; }

    void reserve(Priority bucket, int capacity) { container[bucket].reserve(capacity); }

protected:
    Priority updateLastNonEmptyIndex() {
        for (; lastNonEmptyIndex >= firstNonEmptyIndex; --lastNonEmptyIndex) {
            if (!container[lastNonEmptyIndex].empty()) { return lastNonEmptyIndex; }
        }
        return InvalidIndex;
    }


    // `container[bucket][i]` is the i_th item in the bucket.
    ContainerType container;
    // the index of the first non-empty group.
    Priority firstNonEmptyIndex;
    // the index of the last non-empty group.
    Priority lastNonEmptyIndex;
    // the number of items.
    Priority count;
};


// [RadixBased][SeperatedKey][PositivePriorityOnly]
template<typename T>
class PriorityQueueByBucketL1AutoResizeMinTop : public PriorityQueueByBucketL1FixedSizeMinTop<T> {
public:
    using ContainerType = Vec<Vec<T>>;
    using PriorityQueueByBucketL1FixedSizeMinTop<T>::DefaultHintItemNum;
    using PriorityQueueByBucketL1FixedSizeMinTop<T>::container;
    using PriorityQueueByBucketL1FixedSizeMinTop<T>::PriorityQueueByBucketL1FixedSizeMinTop;


    // put the item in proper place regarding its precedence.
    void push(const T &e, Priority key) {
        if (key >= container.size()) { container.resize(key + 1); } // OPT[szx][1]: more aggressive reallocation?
        PriorityQueueByBucketL1FixedSizeMinTop<T>::push(e, key);
    }

    template<typename RandomNumGenerator>
    void push(const T &e, Priority key, RandomNumGenerator &rgen) {
        if (key >= container.size()) { container.resize(key + 1); }
        PriorityQueueByBucketL1FixedSizeMinTop<T>::push(e, key, rgen);
    }
};

template<typename T>
using PriorityQueueByBucketL1 = PriorityQueueByBucketL1FixedSizeMinTop<T>;


template<typename T, typename LessKeyPred = LessKeyCmp<T>>
using PriorityQueue = PriorityQueueBySet<T, LessKeyPred>;

}
}


#endif // CN_HUST_HAFV_UTIL_PRIORITY_QUEUE_H
