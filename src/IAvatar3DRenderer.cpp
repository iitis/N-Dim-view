#include "IAvatar3DRenderer.h"
#include <qopengl.h>
#include "Avatar3D.h"

void IAvatar3DRenderer::renderSelf(const CBaseObject* _obj)
{
	Avatar3D* obj = (Avatar3D*)_obj;
	obj->draw();
}


