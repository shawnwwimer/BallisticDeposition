#include "FileWriting.h"

#define COMPRESSION_LEVEL (Z_DEFAULT_COMPRESSION)


int writeFileToZipCTS(const char* zipname, const char* filename)
{
	int error = ZIP_ERRNO;

	zipFile zf = zipOpen64(zipname, APPEND_STATUS_ADDINZIP);

	// Check for a valid zipfile
	if (zf == NULL) {
		return ZIP_BADZIPFILE;
	}

	// Attempt to open the file
	std::fstream file(filename, std::ios::binary | std::ios::in);
	if (file.is_open()) {
		// Get size of file
		file.seekg(0, std::ios::end);
		size_t size = file.tellg();
		file.seekg(0, std::ios::beg);

		// Read in file 
		// TODO: might want to chunk it since some files may be large
		std::vector<char> buffer(size);
		if (size == 0 || file.read(&buffer[0], size)) {
			// Initialize the parameters for the local header
			zip_fileinfo zi = { 0 };

			// Open it inside the zip for writing
			if (Z_OK == zipOpenNewFileInZip64(zf, filename, &zi, NULL, 0, NULL, 0, NULL, Z_DEFLATED, COMPRESSION_LEVEL, 1)) {
				// Write in zip
				if (zipWriteInFileInZip(zf, size == 0 ? "" : &buffer[0], size)) {
					error = ZIP_ERRNO;
				}

				// Close it inside zip
				if (zipCloseFileInZip(zf) == ZIP_OK) {
					error = ZIP_OK;
				}
				// Close file
				file.close();
			}
		}
	}

	// Try to close zip
	if (zipClose(zf, NULL)) {
		return ZIP_ERRNO;
	}
	return error;
}