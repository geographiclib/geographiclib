# Script to build static library prefab for android.

## requirement :

- Android NDK r24 : https://developer.android.com/ndk/downloads#r24-downloads

- Google Prefab cli : https://github.com/google/prefab

- Maven Cli

## build and publish :

run build.sh script to publish aar package to local maven repository and ossrh public repository.

to upload to Sonatype OSSRH repository you need to have 'ossrh' server setting in '~/.m2/setting.xml' file.

## Integration

### build.gradle

```gradle
android {
    ...
    buildFeatures {
        ...
        prefab true
    }
}

dependencies {
    ...
    implementation "net.sf.geographiclib:geographiclib-prefab:2.1.0"
}
```

### CMakeLists.txt

```cmake
add_library(app SHARED app.cpp)

# Add these two lines.
find_package (geographiclib REQUIRED CONFIG)
target_link_libraries(app geographiclib::geographiclib)
```
