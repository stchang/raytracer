#lang racket

(require racket/draw
         plot/utils ; vector fns
         racket/generic)

(define WIDTH 640)
(define HEIGHT 480)

;; color pixels -- first number is alpha
(define red-pix (bytes 255 255 0 0))
(define green-pix (bytes 255 0 255 0))
(define blue-pix (bytes 255 0 0 255))
(define black-pix (bytes 255 0 0 0))
(define clear-pix (bytes 0 0 0 0))
(define white-pix (bytes 255 255 255 255))

;; A Point is a vector #(x y z).

;; A Light is a (Light Point Number)
;;  The Point is the position of the light.
;;  The Number is the light's intensity.
(struct Light (pos intensity))

;; A Ray is a Point that represents a unit vector originating from #(0,0,0).

;; A Cam is a Ray.
;; We assume the Cam is at the origin #(0,0,0). The Ray representing the Cam
;;  indicates the Cam's direction.

;; mk-cam : Number Number Number -> Cam
(define (mk-cam x y z) (vnormalize (vector x y z)))

;; these constants are used when calculating rays
;; -- needed to maintain proper aspect ratio
;;    bc both x and yoffset are multiplied by unit vector
;; ie if WIDTH > HEIGHT, then 0 < xoffset < 1
;;                        and 0 < yoffset < HEIGHT/WIDTH
(define WIDTH-SCALE (if (> WIDTH HEIGHT) 1 (/ WIDTH HEIGHT)))
(define HEIGHT-SCALE (if (> HEIGHT WIDTH) 1 (/ HEIGHT WIDTH)))
(define (mk-ray x y cam right up)
  (vnormalize 
   (v+ cam (v+ (v* right (* (- (/ x WIDTH) 0.5) WIDTH-SCALE))
               (v* up (* (- (/ y HEIGHT) 0.5) HEIGHT-SCALE))))))

;(define rays
;  (for*/list ([x (in-range WIDTH)] [y (in-range HEIGHT)])
;    (mk-ray x y)))
;;    (define xoffset (/ x DIV))
;;    (define yoffset (/ y DIV))
;;    (vnormalize (v+ cam (v+ (v* right xoffset)
;;                            (v* up yoffset))))))

(define-generics shape
  ;; intersect : Shape Point [Point] -> Number
  ;; Returns the distance between start and the intersection point of the vector
  ;; represented by start and dir with the given shape, if that vector 
  ;; infinitely extended. Returns +inf.0. If no starting point is specified, 
  ;; then #(0,0,0) is assumed.
  (intersect shape dir [start])
  ;; get-norm : Shape Point -> Ray
  ;; Returns the Ray that is perpendicular to the Shape, at the given point.
  (get-norm shape pt))

(struct Shape (color))

(struct Sphere (center radius) 
  #:super struct:Shape
  #:methods gen:shape
  [(define (intersect s dir [pos #(0 0 0)])
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
  [(define (intersect this-plane dir [pos #(0 0 0)])
     (define pt (Plane-pt this-plane))
     (define norm (Plane-norm this-plane))
     (define denom (vdot dir norm))
     (if (zero? denom) +inf.0
         (let* ([d1 (vdot pos norm)]
                [d2 (vdot pt norm)]
                [dist (/ (- d2 d1) denom)])
           (if (< dist 0) +inf.0 dist))))
  (define (get-norm this-plane pt) (Plane-norm this-plane))])
         

;(define pix (make-bytes 32))
;
;(send scene get-argb-pixels 0 0 2 2 pix)
;
;pix 

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
    (< (intersect other-shape (vnormalize vec-to-light) isectpt)
       dist-to-light)))

(define (render cam lght shapes filename)
  
  (define scene (make-bitmap WIDTH HEIGHT))

  ;; right, with respect to the direction of cam
  (define right (vnormalize (vcross #(0 -1 0) cam)))
  ;; up, with respect to the direction of cam
  (define up (vnormalize (vcross right cam)))

  (for* ([x (in-range WIDTH)]
         [y (in-range HEIGHT)])
    (define ray (mk-ray x y cam right up))
    (define min-dist-shape (argmin (λ (sh) (intersect sh ray)) shapes))
    (define min-dist (intersect min-dist-shape ray)) ;; duplicate calc
    (send 
     scene 
     set-argb-pixels x (- HEIGHT y 1) 1 1
     (if (= min-dist +inf.0)
         clear-pix
         (let* ([isectpt (v* ray min-dist)]
                [normray (get-norm min-dist-shape isectpt)]
                [lightreflectray (vnormalize (v- (Light-pos lght) isectpt))])
           (define reflect-intensity (abs (vdot lightreflectray normray)))
           (cond [(> reflect-intensity 1) clear-pix]
                 [(blocked? isectpt min-dist-shape shapes lght) black-pix]
                 [else
                  (define intensity (* (Light-intensity lght) reflect-intensity))
                  (scale-color (Shape-color min-dist-shape) intensity)])))))

  (send scene save-file filename 'png))

(define cam (mk-cam 0 0 -1))

(define sph1 (Sphere red-pix #(4 12 10) 3))
(define sph2 (Sphere red-pix #(0 8 8) 1))
(define sph3 (Sphere red-pix #(2 2 -10) 1))

(define shapes (list sph3))

(define lght (Light #(0 10 -2) 1))

;(render cam lght shapes "a.png")
(render (mk-cam 5 -5 10) 
        (Light #(15 15 -10) 1)
        (list (Plane white-pix #(5 -6 10) (vnormalize #(0 1 0)))
              (Sphere red-pix #(1 -5 10) 1)
              (Sphere green-pix #(5 -5 10) 1)
              (Sphere blue-pix #(9 -5 10) 1))
        "b.png")