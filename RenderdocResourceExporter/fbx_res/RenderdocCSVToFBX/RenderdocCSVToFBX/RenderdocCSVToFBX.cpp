// RenderdocCSVToFBX.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <set>

#include <fbxsdk.h>

#include "CSVFile.h"
#include "CommonMath.h"

using namespace std;

struct MeshVertex
{
	float x, y, z;

	float nx, ny, nz;
	float tx, ty, tz, tw;

	float u, v;
	float u2, v2;
	float u3, v3;
};

void ConvertCSV2FBX(const char* sourceCSVFile, 
	bool export_normal, bool export_tangent, 
	bool export_uv, bool export_uv2, bool export_uv3)
{
	FbxManager* sdkManager;
	FbxIOSettings* ioSettings;
	FbxScene* scene;

	sdkManager = FbxManager::Create();

	ioSettings = FbxIOSettings::Create(sdkManager, IOSROOT);
	ioSettings->SetBoolProp(EXP_FBX_MATERIAL, true);
	ioSettings->SetBoolProp(EXP_FBX_TEXTURE, true);
	ioSettings->SetBoolProp(EXP_FBX_EMBEDDED, true);
	ioSettings->SetBoolProp(EXP_FBX_ANIMATION, true);
	ioSettings->SetBoolProp(EXP_FBX_SHAPE, true);
	ioSettings->SetIntProp(
		EXP_FBX_EXPORT_FILE_VERSION, FBX_FILE_VERSION_7400);

	sdkManager->SetIOSettings(ioSettings);
	scene = FbxScene::Create(sdkManager, "");

	CCSVFile* pSrcFile = new CCSVFile(sourceCSVFile);

	std::map<int, MeshVertex> verticsMap;
	int iVertexID = 0;
	for (int iRow = 0; iRow < pSrcFile->GetRowNum(); iRow++)
	{
		pSrcFile->GetCellValue("IDX", iRow, iVertexID);

		if (verticsMap.find(iVertexID) == verticsMap.end())
		{
			pSrcFile->GetCellValue("POSITION.x", iRow, verticsMap[iVertexID].x);
			pSrcFile->GetCellValue("POSITION.y", iRow, verticsMap[iVertexID].y);
			pSrcFile->GetCellValue("POSITION.z", iRow, verticsMap[iVertexID].z);

			if(export_normal)
			{
				pSrcFile->GetCellValue("NORMAL.x", iRow, verticsMap[iVertexID].nx);
				pSrcFile->GetCellValue("NORMAL.y", iRow, verticsMap[iVertexID].ny);
				pSrcFile->GetCellValue("NORMAL.z", iRow, verticsMap[iVertexID].nz);
			}
			if(export_tangent)
			{
				pSrcFile->GetCellValue("TANGENT.x", iRow, verticsMap[iVertexID].tx);
				pSrcFile->GetCellValue("TANGENT.y", iRow, verticsMap[iVertexID].ty);
				pSrcFile->GetCellValue("TANGENT.z", iRow, verticsMap[iVertexID].tz);
				pSrcFile->GetCellValue("TANGENT.w", iRow, verticsMap[iVertexID].tw);
			}

			if(export_uv)
			{
				pSrcFile->GetCellValue("TEXCOORD0.x", iRow, verticsMap[iVertexID].u);
				pSrcFile->GetCellValue("TEXCOORD0.y", iRow, verticsMap[iVertexID].v);
			}
			if(export_uv2)
			{
				pSrcFile->GetCellValue("TEXCOORD1.x", iRow, verticsMap[iVertexID].u2);
				pSrcFile->GetCellValue("TEXCOORD1.y", iRow, verticsMap[iVertexID].v2);
			}
			if(export_uv3)
			{
				pSrcFile->GetCellValue("TEXCOORD2.x", iRow, verticsMap[iVertexID].u3);
				pSrcFile->GetCellValue("TEXCOORD2.y", iRow, verticsMap[iVertexID].v3);
			}
		}
	}

	int iTotalVerticsCount = verticsMap.size();

	// 用csv的文件名作为mesh的名字
	std::string meshName = sourceCSVFile;
	replace(meshName.begin(), meshName.end(), '/', '\\');
	meshName = meshName.substr(0, meshName.size() - 4);
	int pos = meshName.find_last_of('\\');
	if (pos != std::string::npos)
	{
		meshName = meshName.substr(pos + 1, meshName.size() - pos - 1);
	}

	FbxNode* rootNode = scene->GetRootNode();
	FbxNode* subshapeNode = FbxNode::Create(rootNode, meshName.c_str());

	// convert fbx mesh file
	FbxMesh* meshFbx = FbxMesh::Create(subshapeNode, subshapeNode->GetName());

	FbxGeometryElementNormal* meshNormal = NULL;
	FbxGeometryElementTangent* meshTangent = NULL;
	FbxGeometryElementUV* meshUV = NULL;
	FbxGeometryElementUV* meshUV2 = NULL;
	FbxGeometryElementUV* meshUV3 = NULL;

	if(export_normal)
	{
		meshNormal = meshFbx->CreateElementNormal();
	}
	if(export_tangent)
	{
		meshTangent = meshFbx->CreateElementTangent();
	}

	if(export_uv)
	{
		meshUV = meshFbx->CreateElementUV("UV");
	}
	if(export_uv2)
	{
		meshUV2 = meshFbx->CreateElementUV("UV1");
	}
	if(export_uv3)
	{
		meshUV3 = meshFbx->CreateElementUV("UV2");
	}

	meshFbx->InitControlPoints(iTotalVerticsCount);

	if(export_normal)
	{
		meshNormal->SetMappingMode(FbxGeometryElementNormal::eByControlPoint);
		meshNormal->SetReferenceMode(FbxGeometryElementNormal::eDirect);
	}
	if(export_tangent)
	{
		meshTangent->SetMappingMode(FbxGeometryElementTangent::eByControlPoint);
		meshTangent->SetReferenceMode(FbxGeometryElementTangent::eDirect);
	}

	if(export_uv)
	{
		meshUV->SetMappingMode(FbxGeometryElementUV::eByControlPoint);
		meshUV->SetReferenceMode(FbxGeometryElementUV::eDirect);
	}
	if(export_uv2)
	{
		meshUV2->SetMappingMode(FbxGeometryElementUV::eByControlPoint);
		meshUV2->SetReferenceMode(FbxGeometryElementUV::eDirect);
	}
	if(export_uv3)
	{
		meshUV3->SetMappingMode(FbxGeometryElementUV::eByControlPoint);
		meshUV3->SetReferenceMode(FbxGeometryElementUV::eDirect);
	}

	FbxVector4* meshVectors = meshFbx->GetControlPoints();
	Matrix44 matRot;
	MatrixRotationZ(&matRot, -FLT_DTOR(0));

	for (int index = 0; index < verticsMap.size(); index++)
	{
		Vector3 Vertex(
			verticsMap[index].x * 100, verticsMap[index].y * 100, verticsMap[index].z * 100);
		Vec3TransformCoord(&Vertex, &Vertex, &matRot);
		meshVectors[index].Set(Vertex.x, Vertex.y, Vertex.z);

		if(export_normal)
		{
			Vector3 VertexNormal(
				verticsMap[index].nx, verticsMap[index].ny, verticsMap[index].nz);
			Vec3TransformNormal(&VertexNormal, &VertexNormal, &matRot);

			meshNormal->GetDirectArray().Add(
				FbxVector4(VertexNormal.x, VertexNormal.y, VertexNormal.z, 0));
		}
		if(export_tangent)
		{
			Vector4 VertexTangent(
				verticsMap[index].tx, verticsMap[index].ty, verticsMap[index].tz, verticsMap[index].tw);
			VertexTangent = VertexTangent * matRot;

			meshTangent->GetDirectArray().Add(
				FbxVector4(VertexTangent.x, VertexTangent.y, VertexTangent.z, VertexTangent.w));
		}
		
		if(export_uv)
		{
			meshUV->GetDirectArray().Add(
				FbxVector2(verticsMap[index].u, verticsMap[index].v));
		}
		if(export_uv2)
		{
			meshUV2->GetDirectArray().Add(
				FbxVector2(verticsMap[index].u2, verticsMap[index].v2));
		}
		if(export_uv3)
		{
			meshUV3->GetDirectArray().Add(
				FbxVector2(verticsMap[index].u3, verticsMap[index].v3));
		}
	}

	int iFaceID = 0;
	for (int iRow = 0; iRow < pSrcFile->GetRowNum(); iRow += 3)
	{
		meshFbx->BeginPolygon(0);

		pSrcFile->GetCellValue("IDX", iRow, iFaceID);
		meshFbx->AddPolygon(iFaceID);
		pSrcFile->GetCellValue("IDX", iRow + 1, iFaceID);
		meshFbx->AddPolygon(iFaceID);
		pSrcFile->GetCellValue("IDX", iRow + 2, iFaceID);
		meshFbx->AddPolygon(iFaceID);

		meshFbx->EndPolygon();
	}

	subshapeNode->SetNodeAttribute(meshFbx);
	subshapeNode->SetShadingMode(FbxNode::eTextureShading);

	// save fbx
	std::string fbxFilePath = sourceCSVFile;
	fbxFilePath = fbxFilePath.substr(0, fbxFilePath.size() - 3);
	fbxFilePath += "fbx";

	int startIndex = fbxFilePath.find_last_of('\\');
	std::string fileName = fbxFilePath.substr(startIndex + 1);

	FbxExporter* exporter = FbxExporter::Create(sdkManager, fileName.c_str());
	if (!exporter->Initialize(fbxFilePath.c_str(), -1, ioSettings))
	{
		fprintf(stderr, "Failed to initialize FBX exporter\n");
		exporter->Destroy();
		return;
	}
	exporter->SetFileExportVersion(FBX_2014_00_COMPATIBLE);

	if (!exporter->Export(scene))
	{
		fprintf(stderr, "Failed to produce FBX file\n");
		exporter->Destroy();
		return;
	}

	exporter->Destroy();
}

int main(int argc, char* argv[])
{
	bool export_normal = false;
	bool export_tangent = false;
	bool export_uv = false;
	bool export_uv2 = false;
	bool export_uv3 = false;

	if (strncmp(argv[2], "1", 1) == 0)
	{
		export_normal = true;
	}
	if (strncmp(argv[3], "1", 1) == 0)
	{
		export_tangent = true;
	}
	if (strncmp(argv[4], "1", 1) == 0)
	{
		export_uv = true;
	}
	if (strncmp(argv[5], "1", 1) == 0)
	{
		export_uv2 = true;
	}
	if (strncmp(argv[6], "1", 1) == 0)
	{
		export_uv3 = true;
	}

	ConvertCSV2FBX(argv[1], 
		export_normal, export_tangent, 
		export_uv, export_uv2, export_uv3);
}
