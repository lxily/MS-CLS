////////////////////////////////
/// 简介:	1.	固定数值范围的整数集合的增删改查操作.
/// 		2.	对于每个整数 i, 总是满足 `(i == item[index[i]])`.
///
/// 备注:	1.	
////////////////////////////////

#ifndef CN_HUST_HAFV_UTIL_CONSECUTIVE_ID_SET_H
#define CN_HUST_HAFV_UTIL_CONSECUTIVE_ID_SET_H


#include <algorithm>

#include "Typedef.h"
#include "Exception.h"
#include "Math.h"


namespace hust {
namespace util {

class ConsecutiveIdSet {
    static constexpr bool SafetyCheck = false;
public:
    using Index = ID; // EXT[szx][4]: make sure no integer overflow! make it a template parameter?
    using Item = Index;


    static constexpr Index InvalidIndex = -1;


    ConsecutiveIdSet(Index capacity, Item minValue = 0)
        : lowerBound(minValue), upperBound(minValue + capacity), itemNum(0),
        items(capacity), index(capacity, InvalidIndex) {
    }


    // if the item has been inserted, return true.
    bool isItemExist(Item e) const { return (index[e - lowerBound] != InvalidIndex); }
    // if the item is in range [min, max], return true.
    bool isItemValid(Item e) const { return math::isInRange(e, lowerBound, upperBound); }
    // if the index is in range [0, itemNum), return true.
    bool isIndexValid(Index i) const { return math::isInRange<Index>(i, 0, itemNum); }

    // return item at the index of i.
    Item itemAt(Index i) const {
        if (SafetyCheck && !isIndexValid(i)) { throw IndexOutOfRangeException(); }
        return items[i];
    }
    // return index of item e.
    Index indexOf(Item e) const {
        if (SafetyCheck && !isItemValid(e)) { throw IndexOutOfRangeException(); }
        if (SafetyCheck && !isItemExist(e)) { throw ItemNotExistException(); }
        return index[e - lowerBound];
    }
    // return the last item.
    Item back() const { return items[size() - 1]; }

    // insert item e and update index.
    void insert(Item e) {
        if (SafetyCheck && !isItemValid(e)) { throw IndexOutOfRangeException(); }
        if (SafetyCheck && isItemExist(e)) { throw DuplicateItemException(); }
        items[itemNum] = e;
        index[e - lowerBound] = itemNum;
        ++itemNum;
    }
    // remove item and update index.
    // never call erase during forward traversing, but it is usually ok for backward traversing.
    void eraseItem(Item e) {
        if (SafetyCheck && !isItemValid(e)) { throw IndexOutOfRangeException(); }
        if (SafetyCheck && !isItemExist(e)) { throw ItemNotExistException(); }
        Index i = index[e - lowerBound];
        --itemNum;
        index[items[itemNum] - lowerBound] = i;
        items[i] = items[itemNum];
        index[e - lowerBound] = InvalidIndex;
    }
    // remove index and update item.
    // never call erase during forward traversing, but it is usually ok for backward traversing.
    void eraseIndex(Index i) {
        if (SafetyCheck && !isIndexValid(i)) { throw IndexOutOfRangeException(); }
        Item e = items[i];
        --itemNum;
        index[items[itemNum] - lowerBound] = i;
        items[i] = items[itemNum];
        index[e - lowerBound] = InvalidIndex;
    }
    // remove the last item.
    Item pop() {
        Item e = items[--itemNum];
        index[e - lowerBound] = InvalidIndex;
        return e;
    }

    // return number of inserted items.
    Index size() const { return itemNum; }

    bool empty() const { return size() <= 0; }

    // keep the capacity but invalidate all items and reset the size.
    void clear(bool sparseData = false) {
        if (sparseData) {
            for (Index i = 0; i < itemNum; ++i) { index[items[i] - lowerBound] = InvalidIndex; }
        } else {
            std::fill(index.begin(), index.end(), InvalidIndex);
        }
        itemNum = 0;
    }

    //const Arr<Item>& getItems() const { return items; } // TODO[szx][0]: handle the invalid items in [itemNum, upperBound).
    const Vec<Index>& getIndices() const { return index; }

protected:
    Item lowerBound; // min value of items.
    Item upperBound; // max value of items.

    Index itemNum; // current number of items.
    Vec<Item> items; // items value, itemNum valid items in it.
    Vec<Index> index; // items index in item.
};

}
}


#endif // CN_HUST_HAFV_UTIL_CONSECUTIVE_ID_SET_H