#pragma once
#include <io.h>
#include <direct.h>
#include <windows.h>
#include <cstring>
#include <string>
#include <fstream>
using namespace std;

//useless??? but work on codeblocks
////判断文件夹是否存在
//bool checkFolderExist(const string &strPath)
//{
//	WIN32_FIND_DATA wfd;
//	bool rValue = false;
//	HANDLE hFind = FindFirstFile(LPCWSTR(strPath.c_str()), &wfd);
//	if ((hFind) != INVALID_HANDLE_VALUE && (wfd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
//		rValue = true;
//	}
//	FindClose(hFind);
//	return rValue;
//}
//
////判断文件是否存在
//bool checkFileExist(const string &filePath){
//	WIN32_FIND_DATA FindFileData;
//	HANDLE hFind;
//
//	hFind = FindFirstFile(LPCWSTR(filePath.c_str()), &FindFileData);
//	if (hFind == INVALID_HANDLE_VALUE) {
//		printf("Invalid File Handle. Get Last Error reports %d.", GetLastError());
//		return false;
//	}
//	else {
//		printf("The first file found is %s.", FindFileData.cFileName);
//		FindClose(hFind);
//		return true;
//	}
//}

//判断文件夹是否存在
bool checkFolderExist(const string &folderPath)
{
	return _access(folderPath.c_str(), 0) != -1;
}

bool makeSureDirectoryPathExists(const string &folderPath) {
	bool exist = checkFolderExist(folderPath);
	if (!exist) system(("md " + folderPath).c_str());
	return true;
}

bool checkFileExist(const string &filePath) {
	fstream _file;
	_file.open(filePath, ios::in);
	return (!!_file);
}

