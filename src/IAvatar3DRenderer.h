#include "../../../dpVision/src/renderers/IObjectRenderer.h"

class CBaseObject;
class CModel3D;

class IAvatar3DRenderer : public IObjectRenderer {
public:
    virtual void renderSelf(const CBaseObject* _obj) override;
private:
    void renderBoundingBox(CModel3D* _obj);
    void renderAxes();
};
