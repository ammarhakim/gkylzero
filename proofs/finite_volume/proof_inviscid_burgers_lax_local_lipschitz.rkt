#lang racket

(require "../prover_core.rkt")

 (prove-lax-friedrichs-scalar-1d-local-lipschitz '#hash((cons-expr . u) (flux-expr . (* 0.5 u u)) (max-speed-expr . (abs u)) (name . "burgers") (parameters . ())) #:cfl 0.95 #:init-func '(cond ((< (abs x) 1.0) 3.0) (else -1.0)) #:nx 200 #:t-final 0.5 #:x0 -3.0 #:x1 3.0)
  (is-real '(cond ((< (abs x) 1.0) 3.0) (else -1.0)) '(u) '())
   (is-real 3.0 '(u) '())
   #t
  (is-real -1.0 '(u) '())
  #t
  (symbolic-diff '(* 0.5 u u) 'u)
   (symbolic-diff 0.5 'u)
   0.0
   (symbolic-diff 'u 'u)
   1.0
   (symbolic-diff 'u 'u)
   1.0
  '(+ (* 0.0 u u) (* 0.5 1.0 u) (* 0.5 u 1.0))
  (symbolic-simp '(+ (* 0.0 u u) (* 0.5 1.0 u) (* 0.5 u 1.0)))
   (symbolic-simp-rule '(+ (* 0.0 u u) (* 0.5 1.0 u) (* 0.5 u 1.0)))
   '(+ (+ (* 0.0 u u) (* 0.5 1.0 u)) (* 0.5 u 1.0))
  (symbolic-simp '(+ (+ (* 0.0 u u) (* 0.5 1.0 u)) (* 0.5 u 1.0)))
   (symbolic-simp-rule '(+ (+ (* 0.0 u u) (* 0.5 1.0 u)) (* 0.5 u 1.0)))
   '(+ (* 0.0 u u) (+ (* 0.5 1.0 u) (* 0.5 u 1.0)))
  (symbolic-simp '(+ (* 0.0 u u) (+ (* 0.5 1.0 u) (* 0.5 u 1.0))))
   (symbolic-simp-rule '(+ (* 0.0 u u) (+ (* 0.5 1.0 u) (* 0.5 u 1.0))))
    (symbolic-simp-rule '(* 0.0 u u))
    '(* (* 0.0 u) u)
    (symbolic-simp-rule '(+ (* 0.5 1.0 u) (* 0.5 u 1.0)))
     (symbolic-simp-rule '(* 0.5 1.0 u))
     '(* (* 0.5 1.0) u)
     (symbolic-simp-rule '(* 0.5 u 1.0))
     '(* (* 0.5 u) 1.0)
    '(+ (* (* 0.5 1.0) u) (* (* 0.5 u) 1.0))
   '(+ (* (* 0.0 u) u) (+ (* (* 0.5 1.0) u) (* (* 0.5 u) 1.0)))
  (symbolic-simp '(+ (* (* 0.0 u) u) (+ (* (* 0.5 1.0) u) (* (* 0.5 u) 1.0))))
   (symbolic-simp-rule '(+ (* (* 0.0 u) u) (+ (* (* 0.5 1.0) u) (* (* 0.5 u) 1.0))))
    (symbolic-simp-rule '(* (* 0.0 u) u))
    '(* 0.0 (* u u))
    (symbolic-simp-rule '(+ (* (* 0.5 1.0) u) (* (* 0.5 u) 1.0)))
     (symbolic-simp-rule '(* (* 0.5 1.0) u))
     '(* 0.5 (* 1.0 u))
     (symbolic-simp-rule '(* (* 0.5 u) 1.0))
     '(* 0.5 (* u 1.0))
    '(+ (* 0.5 (* 1.0 u)) (* 0.5 (* u 1.0)))
   '(+ (* 0.0 (* u u)) (+ (* 0.5 (* 1.0 u)) (* 0.5 (* u 1.0))))
  (symbolic-simp '(+ (* 0.0 (* u u)) (+ (* 0.5 (* 1.0 u)) (* 0.5 (* u 1.0)))))
   (symbolic-simp-rule '(+ (* 0.0 (* u u)) (+ (* 0.5 (* 1.0 u)) (* 0.5 (* u 1.0)))))
    (symbolic-simp-rule '(* 0.0 (* u u)))
    0.0
    (symbolic-simp-rule '(+ (* 0.5 (* 1.0 u)) (* 0.5 (* u 1.0))))
     (symbolic-simp-rule '(* 0.5 (* 1.0 u)))
     '(* 0.5 u)
     (symbolic-simp-rule '(* 0.5 (* u 1.0)))
      (symbolic-simp-rule 0.5)
      0.5
      (symbolic-simp-rule '(* u 1.0))
      '(* 1.0 u)
     '(* 0.5 (* 1.0 u))
    '(+ (* 0.5 u) (* 0.5 (* 1.0 u)))
   '(+ 0.0 (+ (* 0.5 u) (* 0.5 (* 1.0 u))))
  (symbolic-simp '(+ 0.0 (+ (* 0.5 u) (* 0.5 (* 1.0 u)))))
   (symbolic-simp-rule '(+ 0.0 (+ (* 0.5 u) (* 0.5 (* 1.0 u)))))
   '(+ (* 0.5 u) (* 0.5 (* 1.0 u)))
  (symbolic-simp '(+ (* 0.5 u) (* 0.5 (* 1.0 u))))
   (symbolic-simp-rule '(+ (* 0.5 u) (* 0.5 (* 1.0 u))))
    (symbolic-simp-rule '(* 0.5 u))
     (symbolic-simp-rule 0.5)
     0.5
     (symbolic-simp-rule 'u)
     'u
    '(* 0.5 u)
    (symbolic-simp-rule '(* 0.5 (* 1.0 u)))
    '(* 0.5 u)
   '(+ (* 0.5 u) (* 0.5 u))
  (symbolic-simp '(+ (* 0.5 u) (* 0.5 u)))
   (symbolic-simp-rule '(+ (* 0.5 u) (* 0.5 u)))
   '(* (+ 0.5 0.5) u)
  (symbolic-simp '(* (+ 0.5 0.5) u))
   (symbolic-simp-rule '(* (+ 0.5 0.5) u))
    (symbolic-simp-rule '(+ 0.5 0.5))
    1.0
    (symbolic-simp-rule 'u)
    'u
   '(* 1.0 u)
  (symbolic-simp '(* 1.0 u))
   (symbolic-simp-rule '(* 1.0 u))
   'u
  (symbolic-simp 'u)
   (symbolic-simp-rule 'u)
   'u
  'u
  (symbolic-diff 'u 'u)
  1.0
  (symbolic-simp 1.0)
   (symbolic-simp-rule 1.0)
   1.0
  1.0
  (is-non-negative 1.0 '())
  #t
 #t
