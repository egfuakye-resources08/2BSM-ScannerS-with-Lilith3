#include "ScannerS/DataMap.hpp"

#include <utility>

namespace ScannerS {

const DataMap::mapped_type &DataMap::operator[](const key_type &key) const {
  auto res = data_.find(key);
  if (res == data_.end())
    throw DataMapError("Unknown key " + key);
  return res->second;
}

void DataMap::Store(const key_type &key, double value) {
  auto res = data_.try_emplace(key, value);
  if (!res.second)
    throw DataMapError("Can't store, key " + key + " already exists.");
}

void DataMap::Store(key_type &&key, double value) {
  auto res = data_.try_emplace(std::move(key), value);
  if (!res.second)
    throw DataMapError("Can't store, key " + key + " already exists.");
}

void DataMap::Merge(Map &&source) {
  data_.merge(source);
  for (auto [key, value] : source) {
    throw DataMapError("Entry {" + key + ", " + std::to_string(value) +
                       "} can't be merged. Exists with value " +
                       std::to_string(data_.at(key)));
  }
}

} // namespace ScannerS
