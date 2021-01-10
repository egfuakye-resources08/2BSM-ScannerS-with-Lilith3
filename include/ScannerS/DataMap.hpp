#pragma once

#include <map>
#include <stdexcept>
#include <string>

namespace ScannerS {

//! An associative container where elemnts cannot be modified after insertion
class DataMap {
public:
  using Map = std::map<std::string, double>;

  using key_type = Map::key_type;       //!< type of the keys
  using mapped_type = Map::mapped_type; //!< type of the mapped values

  /** Return the value corresponding to the key.
   * If key does not exist this throws a DataMapError.
   *
   * @param key identifies the Element
   * @return const mapped_type& the value corresponding to key
   */
  const mapped_type &operator[](const key_type &key) const;

  /** Adds a new element {key : value}.
   * If key already exists this throws a DataMapError.
   *
   * @param key identifies the Element
   * @param value the value corresponding to key
   */
  void Store(const key_type &key, mapped_type value);
  //! Overload for Rvalue keys @copydoc Store()
  void Store(key_type &&key, mapped_type value);

  /** Stores all entries of the source.
   * Throws a DataMapError it one of the keys already exists.
   *
   * @param source
   */
  void Merge(Map &&source);

  //! Const begin iterator
  auto begin() const noexcept { return data_.begin(); }
  //! Const end iterator
  auto end() const noexcept { return data_.end(); }

private:
  Map data_ = Map{};
};

//! Error class thrown by DataMap.
class DataMapError : public std::runtime_error {
public:
  //! Constructor taking a message
  DataMapError(const std::string &message) : runtime_error{message} {}
};

} // namespace ScannerS
