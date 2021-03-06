#ifndef CMU462_DYNAMICSCENE_MESHVIEW_H
#define CMU462_DYNAMICSCENE_MESHVIEW_H

namespace CMU462 { namespace DynamicScene {

/*
  An interface used to give access to mesh-specific functionality (edge
  flipping, subdivision, etc) without depending on a specific implementation.
*/
class MeshView {
 public:
  /*
    Operations to do with a specific feature on the mesh.
  */
  virtual void collapse_selected_edge() = 0;
  virtual void flip_selected_edge() = 0;
  virtual void split_selected_edge() = 0;
  /*
    Operations to do with the mesh as a whole.
  */
  virtual void upsample() = 0;
  virtual void downsample() = 0;
  virtual void resample() = 0;
  virtual void drag_selection_normal(float dx, float dy,
                                     const Matrix4x4& worldTo3DH) = 0;
  virtual void init_animation(int type) = 0;
  virtual void animate(int type) = 0;
  virtual void propogate() = 0;

};

} // namespace DynamicScene
} // namespace CMU462

#endif //CMU462_DYNAMICSCENE_MESHVIEW_H
