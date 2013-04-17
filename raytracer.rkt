#lang racket

(require racket/draw
         plot/utils ; vector fns
         racket/generic)

;; color pixels -- first number is alpha
(define red-pixel (bytes 255 255 0 0))
(define green-pixel (bytes 255 0 255 0))
(define blue-pixel (bytes 255 0 0 255))
(define black-pixel (bytes 255 0 0 0))
(define clear-pixel (bytes 0 0 0 0))
(define white-pixel (bytes 255 255 255 255))


;; A Point is a (Racket) vector #(x y z) representing coordinate (x,y,z).

;; A Dir is a Point representing a direction.
;; The Point represents a unit length (Euclidean) vector with origin (0,0,0).

;; A Light is a (Light pos intensity) : Point Number -> Light
;; representing a light source where
;; - pos represents the position of the light
;; - intensity represent's the light's intensity.
(struct Light (pos intensity))

;; A Ray is a (Ray pos dir) : Point Dir -> Ray
;; representing an infinite (Euclidean) vector where
;; - pos is the tail of the Ray
;; - dir is the direction of the Ray.
(struct Ray (pos dir))

;; A Cam is a (Cam pos dir right up) : Point Dir Dir Dir -> Cam
;; representing a camera where
;; - pos is the position
;; - dir is the direction where the camera is pointed
;; - right is the right direction relative to the camera's forward direction
;; - up is the camera's relative up direction.
(struct Cam (pos dir right up))

;; mk-cam : Number Number Number -> Cam
;; Takes two points, one for position, and one for direction, and returns a Cam
(define (mk-cam pos look) 
  (define dir (vnormalize (v- look pos)))
  (define right (vnormalize (vcross #(0 -1 0) dir)))
  (define up (vnormalize (vcross right dir)))
  (Cam pos dir right up))


(define-generics shape
  ;; intersect : Shape Ray -> Number
  ;; Returns the distance from the start of the given ray to the shape
  ;; Returns +inf.0 if no intersection.
  (intersect shape ray)
  ;; get-norm : Shape Point -> Ray
  ;; Returns the Ray that is perpendicular to the Shape, at the given point.
  (get-norm shape pt))

(struct Shape (color))

(struct Sphere (center radius) 
  #:super struct:Shape
  #:methods gen:shape
  [(define (intersect s ray)
     (match-define (Ray pos dir) ray)
     (define center (Sphere-center s))
     (define radius (Sphere-radius s))
     (define vec-to-center (v- center pos))
     ;; proj is the projection of center onto ray (forming rt triangle)
     (define proj (vdot vec-to-center dir))
     ;; neg proj means angle between vecs > 90, ie going away from each other
     (if (< proj 0) +inf.0
                ;; r^2 + proj^2 = hypotenuse^2 of a rt triangle
         (let* ([hyp^2 (+ (* radius radius)
                          (* proj proj))]
                ;; (vdot cent cent) = (dist from 0 to center)^2
                [tmp (- hyp^2 (vdot vec-to-center vec-to-center))])
           ;; if this hyp^2 > (dist from 0 to center)^2, then the ray intersects
           ;; ow ray is outside the circle
           (if (< tmp 0) +inf.0
               ;; writing out both sets of pythag eq and subtracting:
               ;; hyp^2 - (vdot cen cen) = tmp = r^2 - opp^2 = adj^2
               ;; where adj and opp are sides of a rt triangle where the hyp is
               ;; the center and the ray intersect point
               ;; so (sqrt tmp) gives adj, which when subtracted from proj gives
               ;; the dist to the intersection point
               ;; (it's easier to see if you draw it out on paper)
               (- proj (sqrt tmp))))))
   (define (get-norm s pt)
     (vnormalize (v- pt (Sphere-center s))))])
               

(struct Plane (pt norm)
  #:super struct:Shape
  #:methods gen:shape
  [(define (intersect this-plane ray)
     (match-define (Ray pos dir) ray)
     (define pt (Plane-pt this-plane))
     (define norm (Plane-norm this-plane))
     (define denom (vdot dir norm))
     (if (zero? denom) +inf.0
         (let* ([d1 (vdot pos norm)]
                [d2 (vdot pt norm)]
                [dist (/ (- d2 d1) denom)])
           (if (< dist 0) +inf.0 dist))))
  (define (get-norm this-plane pt) (Plane-norm this-plane))])
         
;; bs is (bytes alpha r g b) representing a color
;; scale is number in [0,1]
(define (scale-color bs scale)
  (define argb (bytes->list bs))
  (apply bytes 
         (first argb)
         (map 
          (compose inexact->exact round (λ (x) (* scale x)))
          (cdr argb))))

;; returns true if there is a shape between isectpt on given shape and light
(define (blocked? isectpt shape other-shapes lght)
  (define vec-to-light (v- (Light-pos lght) isectpt))
  (define dist-to-light (vmag vec-to-light))
  (for/or ([other-shape other-shapes]
           #:unless (eq? other-shape shape))
    (< (intersect other-shape (Ray isectpt (vnormalize vec-to-light)))
       dist-to-light)))

                               
;; the 0.5 adjustment is so the camera is centered, as opposed to in the
;; top-left of the picture
(define (mk-ray x width width-scale y height height-scale cam)
  (match-define (Cam pos dir right up) cam)
  (Ray
   pos
   (vnormalize 
    (v+ dir (v+ (v* right (* (- (/ x width) 0.5) width-scale))
                (v* up (* (- (/ y height) 0.5) height-scale)))))))

(define (render width height background-color cam lght shapes filename)
  (define scene (make-bitmap width height))
  ;; these constants are used when calculating rays
  ;; -- needed to maintain proper aspect ratio
  (define width-scale (if (> width height) 1 (/ width height)))
  (define height-scale (if (> height width) 1 (/ height width)))
  (for* ([x (in-range width)]
         [y (in-range height)])
    (define ray (mk-ray x width width-scale y height height-scale cam))
    (match-define (Ray ray-pos ray-dir) ray)
    (define min-dist-shape (argmin (λ (sh) (intersect sh ray)) shapes))
    (define min-dist (intersect min-dist-shape ray)) ;; duplicate calc
    (send 
     scene 
     ;; subtract from height so 0 height is bottom instead of top
     set-argb-pixels x (- height y 1) 1 1
     (if (= min-dist +inf.0)
         clear-pixel
         (let* ([isectpt (v+ ray-pos (v* ray-dir min-dist))]
                [normray (get-norm min-dist-shape isectpt)]
                [lightreflectray (vnormalize (v- (Light-pos lght) isectpt))])
           (define reflect-intensity (abs (vdot lightreflectray normray)))
           (cond [(> reflect-intensity 1) clear-pixel]
                 [(blocked? isectpt min-dist-shape shapes lght) black-pixel]
                 [else
                  (define intensity (* (Light-intensity lght) reflect-intensity))
                  (scale-color (Shape-color min-dist-shape) intensity)])))))
  (send scene save-file filename 'png))

;(define cam (mk-cam #(0 0 0) #(0 0 -1)))
;(define sph1 (Sphere red-pixel #(4 12 10) 3))
;(define sph2 (Sphere red-pixel #(0 8 8) 1))
;(define sph3 (Sphere red-pixel #(2 2 -10) 1))
;(define shapes (list sph3))
;(define lght (Light #(0 10 -2) 1))
;(render 640 480 clear-pixel cam lght shapes "a.png")

(render 640 480
        clear-pixel
        (mk-cam #(0 0 0) #(5 -5 10))
        (Light #(15 15 -10) 1)
        (list (Plane white-pixel #(5 -6 10) (vnormalize #(0 1 0)))
              (Sphere red-pixel #(1 -5 10) 1)
              (Sphere green-pixel #(5 -5 10) 1)
              (Sphere blue-pixel #(9 -5 10) 1))
        "b.png")