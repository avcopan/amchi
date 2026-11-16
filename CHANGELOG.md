# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- `AmchiError`, the base of the errors that amchi raises, and `UnsupportedElementError` for an element without data (amchi supports the main-group elements from H to Xe).
- The input models `Connectivity` (the bonds, and the breaking and forming bonds of a TS) and `Candidates` (each stereocenter candidate with its stereo-determining neighbors), the `GeometryLike` protocol, and `IllDefinedParityError`.
