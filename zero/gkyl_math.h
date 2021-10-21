#pragma once

// Three component vector
struct  gkyl_vec3 { double x[3]; };

// new vector with zeros 
static inline struct gkyl_vec3
gkyl_vec3_zeros()
{
  return (struct gkyl_vec3) { .x = { 0.0, 0.0, 0.0} };
}

// a+b
static inline struct gkyl_vec3
gkyl_vec3_add(struct gkyl_vec3 a, struct gkyl_vec3 b)
{
  return (struct gkyl_vec3) { .x = { a.x[0]+b.x[0], a.x[1]+b.x[1], a.x[2]+b.x[2] } };
}

// a-b
static inline struct gkyl_vec3
gkyl_vec3_sub(struct gkyl_vec3 a, struct gkyl_vec3 b)
{
  return (struct gkyl_vec3) { .x = { a.x[0]-b.x[0], a.x[1]-b.x[1], a.x[2]-b.x[2] } };
}

// |a|
static inline double
gkyl_vec3_len(struct gkyl_vec3 a)
{
  return sqrt( a.x[0]*a.x[0] + a.x[1]*a.x[1] + a.x[2]*a.x[2] );
}

// normalize a
static inline struct gkyl_vec3
gkyl_vec3_norm(struct gkyl_vec3 a)
{
  double len = gkyl_vec3_len(a);
  return (struct gkyl_vec3) { .x = { a.x[0]/len, a.x[1]/len, a.x[2]/len } };
}

// a \dot b
static inline double
gkyl_vec3_dot(struct gkyl_vec3 a, struct gkyl_vec3 b)
{
  return a.x[0]*b.x[0] + a.x[1]*b.x[1] + a.x[2]*b.x[2];
}

// a \times b
static inline struct gkyl_vec3
gkyl_vec3_cross(struct gkyl_vec3 a, struct gkyl_vec3 b)
{
  return (struct gkyl_vec3) { .x = {
      a.x[1]*b.x[2]-a.x[2]*b.x[1],
      a.x[2]*b.x[0]-a.x[0]*b.x[2],
      a.x[0]*b.x[1]-a.x[1]*b.x[0]
    }
  };
}
