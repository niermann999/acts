
# JSON plugin

## General

The `JSON` plugin facilitates the use of the [nlohmann json](https://github.com/nlohmann/json) library in Acts. It is made available by the `-DACTS_BUILD_PLUGIN_JSON=ON` flag during build configuration and provides the serialization and deserialization of several core Acts data object states, such as material, surfaces or binning data. It is used e.g. during the material mapping step or in the digitization configuration.

In order to automatically read and write a custom object states in json, the nlohmann library needs to be provided with two functions 
```cpp 
to_json(json& j, const object_t& obj) {...}
```
and 
```cpp
from_json(const json& j, object_t& obj) {...}
```
which will be called when the corresponding type is read in or written to a json file. The `JSON` plugin provides a number of implementations that are needed to write out an entire tracking geometry including material (see below).

For the complete documentation of nlohmann json, please see [doc](https://github.com/nlohmann/json).

## Available JSON io implementations

Among the Acts objects that are convertible to/from json via this plugin are:
- {class}`Acts::TrackingGeometry`: In order to write a complete tracking geometry to json and back, several other core classes are made available for json io: The {class}`Acts::TrackingVolume` and the {class}`Acts::Surface`, including all of the {class}`Acts::SurfaceBounds` implementations and its placement via {type}`Acts::Transform3`.
- {class}`Acts::GeometryHierarchyMap`:
- {class}`Acts::BinningData`/{class}`Acts::BinUtility`: 
- {class}`Acts::Material`:

## Usage

In order to read / write a data object from / to a `.json` file, the nlohmann library provides the following syntax:


To then e.g. write a `TrackingVolume` to file, it can 
