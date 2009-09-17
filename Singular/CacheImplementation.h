#ifndef CacheImplementationIncluded
#define CacheImplementationIncluded

template<class KeyClass, class ValueClass>
Cache<KeyClass, ValueClass>::Cache (const int maxEntries, const int maxWeight) {
    _maxEntries = maxEntries;
    _maxWeight = maxWeight;
    _rank.clear();
    _key.clear();
    _value.clear();
    _itKey = _key.end(); // referring to past-the-end element in the list
    _itValue = _value.end(); // referring to past-the-end element in the list
    _weight = 0;
}

template<class KeyClass, class ValueClass>
int Cache<KeyClass, ValueClass>::getWeight() const {
    return _weight;
}

template<class KeyClass, class ValueClass>
int Cache<KeyClass, ValueClass>::getNumberOfEntries() const {
    return _rank.size();
}

template<class KeyClass, class ValueClass>
void Cache<KeyClass, ValueClass>::clear() {
    _rank.clear();
    _key.clear();
    _value.clear();
}

template<class KeyClass, class ValueClass>
Cache<KeyClass, ValueClass>::~Cache() {
    _rank.clear();
    _key.clear();
    _value.clear();
}

template<class KeyClass, class ValueClass>
bool Cache<KeyClass, ValueClass>::hasKey (const KeyClass& key) const {
    _itKey = _key.end(); // referring to past-the-end element in the list
    typename list<KeyClass>::const_iterator itKey;
    _itValue = _value.begin();
    // As _key is a sorted list, the following could actually be implemented
    // in logarithmic time, by bisection. For lists this does not work.
    // But often, we can still terminate the linear loop.
    for (itKey = _key.begin(); itKey != _key.end(); itKey++) {
        int c = key.compare(*itKey);
        if (c == 0) {
            _itKey = itKey;
            return true;
        }
        if (c == -1) return false;
        _itValue++;
    }
    return false;
}

template<class KeyClass, class ValueClass>
ValueClass Cache<KeyClass, ValueClass>::getValue (const KeyClass& key) const {
    if (_itKey == _key.end())
        // _itKey refers to past-the-end element in the list
        // getValue has been called although hasKey produced no match
        assert(false);

    return *_itValue;
}

template<class KeyClass, class ValueClass>
bool Cache<KeyClass, ValueClass>::shrink(const KeyClass& key) {
    // We need to return true iff the pair with given key needed to be erased
    // during the shrinking procedure. So far, we assume no:
    bool result = false;
    // Shrink until both bounds will be met again:
    while (int(_key.size()) > _maxEntries || _weight > _maxWeight) {
        if (deleteLast(key)) result = true;
    }
    return result;
}

template<class KeyClass, class ValueClass>
int Cache<KeyClass, ValueClass>::getMaxNumberOfEntries() const {
    return _maxEntries;
}

template<class KeyClass, class ValueClass>
int Cache<KeyClass, ValueClass>::getMaxWeight() const {
    return _maxWeight;
}

template<class KeyClass, class ValueClass>
bool Cache<KeyClass, ValueClass>::deleteLast(const KeyClass& key) {
    if (_rank.size() == 0) {
        return false; // nothing to do
    };
    // We need to do the next empty loop in order to obtain a forward-iterator pointing to the last
    // entry of _rank.
    // Note: We cannot use rbegin() because we need the iterator for erasing the last entry which is
    // only implemented for forward iterators by std::list.
    list<int>::iterator itRank;
    for (itRank = _rank.begin(); itRank != _rank.end(); itRank++) { }
    itRank--;  // Now, this forward iterator points to the last list entry.
    int deleteIndex = *itRank;   // index of (_key, _value)-pair with worst, i.e., highest _rank
    bool result = false;

    // now delete entries in _key and _value with index deleteIndex
    int k = 0;
    typename list<KeyClass>::iterator itKey;
    typename list<ValueClass>::iterator itValue = _value.begin();
    for (itKey = _key.begin(); itKey != _key.end(); itKey++) {
        if (k == deleteIndex) {
            result = (key.compare(*itKey) == 0);
            break;
        }
        itValue++;
        k++;
    }
    _key.erase(itKey);
    int deleteWeight = itValue->getWeight();
    _value.erase(itValue);

    // adjust total weight of this cache
    _weight -= deleteWeight;

    // now delete last entry of _rank and decrement all those indices
    // in _rank by 1 which were larger than deleteIndex
    _rank.erase(itRank);
    for (itRank = _rank.begin(); itRank != _rank.end(); itRank++) {
        if (*itRank > deleteIndex) *itRank -= 1;
    }

    return result;
}

template<class KeyClass, class ValueClass>
bool Cache<KeyClass, ValueClass>::put (const KeyClass& key, const ValueClass& value) {
    bool keyWasContained = false;
    int oldIndexInKey = -1;
    int newIndexInKey = _key.size();  // default to enter new (key, value)-pair is at the end of the two lists
                                      // only used in the case keyWasContained == false
    int k = 0;
    typename list<KeyClass>::iterator itKey;
    typename list<ValueClass>::iterator itOldValue = _value.begin(); // will later only be used in the case keyWasContained == true
    for (itKey = _key.begin(); itKey != _key.end(); itKey++) {
        int c = key.compare(*itKey);
        if (c == -1) {
            newIndexInKey = k;
            break;
        }
        if (c == 0) {
            keyWasContained = true;
            oldIndexInKey = k;
            break;
        }
        itOldValue++;
        k++;
    }
    int utility = value.getUtility();
    k = 0;
    typename list<ValueClass>::iterator itValue = _value.begin();
    for (itValue = _value.begin(); itValue != _value.end(); itValue++) {
        if (itValue->getUtility() > utility) k++;
    }
    int newIndexInRank = k;

    if (keyWasContained) {
        // There was already a pair of the form (key --> *).

        // adjusting the weight of the cache
        ValueClass oldValue = *itOldValue;
        _weight += value.getWeight() - oldValue.getWeight();

        // overriding old value by argument value
        itOldValue = _value.erase(itOldValue);
        ValueClass myValueCopy = value;
        _value.insert(itOldValue, myValueCopy);

        int oldIndexInRank = -1;
        // oldIndexInRank is to be the position in _rank such that
        // _rank[oldIndexInRank] == oldIndexInKey, i.e.
        // _key[_rank[oldIndexInRank]] == key:
        list<int>::iterator itRank;
        k = 0;
        for (itRank = _rank.begin(); itRank != _rank.end(); itRank++) {
            if (*itRank == oldIndexInKey) {
                oldIndexInRank = k;
            }
            k++;
        }
        // Although the key stays the same, the ranking of the (key --> value) pair can be
        // completely different from before. Thus, we need to repair the entries of _rank:
        if (oldIndexInRank < newIndexInRank) {  // first insert, then erase
            k = 0;
            // insert 'oldIndexInKey' at new position 'newIndexInRank':
            for (itRank = _rank.begin(); itRank != _rank.end(); itRank++) {
                if (k == newIndexInRank) break;
                k++;
            }
            _rank.insert(itRank, oldIndexInKey); // note that this may also insert at position itRank == _rank.end(), i.e.
                                                 // when above loop did not terminate due to 'break' statement
            k = 0;
            // erase 'oldIndexInKey' at old position 'oldIndexInRank':
            for (itRank = _rank.begin(); itRank != _rank.end(); itRank++) {
                if (k == oldIndexInRank) {
                    _rank.erase(itRank);
                    break;
                }
                k++;
            }
        }
        else {  // oldIndexInRank >= newIndexInRank
            if (oldIndexInRank > newIndexInRank) { // first erase, then insert
                k = 0;
                // erase 'oldIndexInKey' at old position 'oldIndexInRank':
                for (itRank = _rank.begin(); itRank != _rank.end(); itRank++) {
                    if (k == oldIndexInRank) {
                        _rank.erase(itRank);
                        break;
                    }
                    k++;
                }
                k = 0;
                // insert 'oldIndexInKey' at new position 'newIndexInRank':
                for (itRank = _rank.begin(); itRank != _rank.end(); itRank++) {
                    if (k == newIndexInRank) {
                        _rank.insert(itRank, oldIndexInKey);
                        break;
                    }
                    k++;
                }
            }
        }
    }
    else {
        // There is no pair of the form (key --> *). We are about to insert a completely new
        // (key, value)-pair.
        // After this branch, we shall have _key[newIndexInKey] = key; _value[newIndexInKey] = value.
        // Note that, after the above computation, newIndexInKey contains the correct target position.
        // Let's make room for the assignment _rank[newIndexInRank] := newIndexInKey:
        list<int>::iterator itRank;
        for (itRank = _rank.begin(); itRank != _rank.end(); itRank++) {
            if (newIndexInKey <= *itRank) {
                *itRank += 1;
            }
        }
        k = 0;
        for (itRank = _rank.begin(); itRank != _rank.end(); itRank++) {
            if (k == newIndexInRank) break;
            k++;
        }
        _rank.insert(itRank, newIndexInKey);
        // let's insert new key and new value at index newIndexInKey:
        itValue = _value.begin();
        k = 0;
        for (itKey = _key.begin(); itKey != _key.end(); itKey++) {
            if (k == newIndexInKey) break;
            itValue++;
            k++;
        }
        KeyClass myKeyCopy = key;
        ValueClass myValueCopy = value;
        _key.insert(itKey, myKeyCopy);
        _value.insert(itValue, myValueCopy);
        // adjusting the total weight of the cache:
        _weight += value.getWeight();
    };
    // We may now have to shrink the cache:
    bool result = shrink(key);  // true iff shrinking deletes the new (key, value)-pair

    assert(_rank.size() == _key.size());
    assert(_rank.size() == _value.size());
    return !result; // true iff the new (key --> value) pair is actually in the cache now
}

template<class KeyClass, class ValueClass>
string Cache<KeyClass, ValueClass>::toString() const {
    char h[10];
    string s = "Cache:";
    s += "\n   entries: ";
    sprintf(h, "%d", getNumberOfEntries()); s += h;
    s += " of at most ";
    sprintf(h, "%d", getMaxNumberOfEntries()); s += h;
    s += "\n   weight: ";
    sprintf(h, "%d", getWeight()); s += h;
    s += " of at most ";
    sprintf(h, "%d", getMaxWeight()); s += h;
    if (_key.size() == 0) {
        s += "\n   no pairs, i.e. cache is empty";
    }
    else {
        int k = 1;
        s += "\n   (key --> value) pairs in ascending order of keys:";
        typename list<KeyClass>::const_iterator itKey;
        typename list<ValueClass>::const_iterator itValue = _value.begin();
        for (itKey = _key.begin(); itKey != _key.end(); itKey++) {
            s += "\n      ";
            sprintf(h, "%d", k); s += h;
            s += ". ";
            s += itKey->toString();
            s += " --> ";
            s += itValue->toString();
            itValue++;
            k++;
        }
        s += "\n   (key --> value) pairs in descending order of ranks:";
        list<int>::const_iterator itRank;
        int r = 1;
        for (itRank = _rank.begin(); itRank != _rank.end(); itRank++) {
           int index = *itRank;
           itValue = _value.begin();
           k = 0;
           for (itKey = _key.begin(); itKey != _key.end(); itKey++) {
               if (k == index) break;
               k++;
               itValue++;
           }
           s += "\n      ";
           sprintf(h, "%d", r); s += h;
           s += ". ";
           s += itKey->toString();
           s += " --> ";
           s += itValue->toString();
           r++;
        }
    }
    return s;
}

template<class KeyClass, class ValueClass>
void Cache<KeyClass, ValueClass>::print() const {
    cout << this->toString();
}

template<class KeyClass, class ValueClass>
Cache<KeyClass, ValueClass>::Cache() { }

template<class KeyClass, class ValueClass>
Cache<KeyClass, ValueClass>::Cache(const Cache& c) {
    _rank = c._rank;
    _value = c._value;
    _key = c._key;
    _weight = c._weight;
    _maxEntries = c._maxEntries;
    _maxWeight = c._maxWeight;
}

#endif
