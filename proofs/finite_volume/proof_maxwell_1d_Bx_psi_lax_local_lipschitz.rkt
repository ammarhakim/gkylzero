#lang racket

(require "../prover_core.rkt")
(require "../prover_vector.rkt")

 (prove-lax-friedrichs-vector2-1d-local-lipschitz '#hash((cons-exprs . (Bx psi)) (flux-exprs . ((* b_fact psi) (* b_fact (* (* c c) Bx)))) (max-speed-exprs . ((abs (* b_fact c)) (abs (* b_fact c)))) (name . "maxwell_Bxpsi") (parameters . ((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))) #:cfl 0.95 #:init-funcs '(0.0 (cond ((< x 0.0) 0.5) (else -0.5))) #:nx 200 #:t-final 1.0 #:x0 -1.5 #:x1 1.5)
  (symbolic-hessian '(* b_fact psi) '(Bx psi))
   (symbolic-gradient '(* b_fact psi) '(Bx psi))
    (symbolic-diff '(* b_fact psi) 'Bx)
     (symbolic-diff 'b_fact 'Bx)
     0.0
     (symbolic-diff 'psi 'Bx)
     0.0
    '(+ (* 0.0 psi) (* b_fact 0.0))
    (symbolic-simp '(+ (* 0.0 psi) (* b_fact 0.0)))
     (symbolic-simp-rule '(+ (* 0.0 psi) (* b_fact 0.0)))
      (symbolic-simp-rule '(* 0.0 psi))
      0.0
      (symbolic-simp-rule '(* b_fact 0.0))
      '(* 0.0 b_fact)
     '(+ 0.0 (* 0.0 b_fact))
    (symbolic-simp '(+ 0.0 (* 0.0 b_fact)))
     (symbolic-simp-rule '(+ 0.0 (* 0.0 b_fact)))
     '(* 0.0 b_fact)
    (symbolic-simp '(* 0.0 b_fact))
     (symbolic-simp-rule '(* 0.0 b_fact))
     0.0
    (symbolic-simp 0.0)
     (symbolic-simp-rule 0.0)
     0.0
    0.0
    (symbolic-diff '(* b_fact psi) 'psi)
     (symbolic-diff 'b_fact 'psi)
     0.0
     (symbolic-diff 'psi 'psi)
     1.0
    '(+ (* 0.0 psi) (* b_fact 1.0))
    (symbolic-simp '(+ (* 0.0 psi) (* b_fact 1.0)))
     (symbolic-simp-rule '(+ (* 0.0 psi) (* b_fact 1.0)))
      (symbolic-simp-rule '(* 0.0 psi))
      0.0
      (symbolic-simp-rule '(* b_fact 1.0))
      '(* 1.0 b_fact)
     '(+ 0.0 (* 1.0 b_fact))
    (symbolic-simp '(+ 0.0 (* 1.0 b_fact)))
     (symbolic-simp-rule '(+ 0.0 (* 1.0 b_fact)))
     '(* 1.0 b_fact)
    (symbolic-simp '(* 1.0 b_fact))
     (symbolic-simp-rule '(* 1.0 b_fact))
     'b_fact
    (symbolic-simp 'b_fact)
     (symbolic-simp-rule 'b_fact)
     'b_fact
    'b_fact
   '(0.0 b_fact)
  (symbolic-jacobian '(0.0 b_fact) '(Bx psi))
   (symbolic-diff 0.0 'Bx)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff 0.0 'psi)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff 'b_fact 'Bx)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff 'b_fact 'psi)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
  '((0.0 0.0) (0.0 0.0))
  (symbolic-hessian '(* b_fact (* (* c c) Bx)) '(Bx psi))
   (symbolic-gradient '(* b_fact (* (* c c) Bx)) '(Bx psi))
    (symbolic-diff '(* b_fact (* (* c c) Bx)) 'Bx)
     (symbolic-diff 'b_fact 'Bx)
     0.0
     (symbolic-diff '(* (* c c) Bx) 'Bx)
      (symbolic-diff '(* c c) 'Bx)
       (symbolic-diff 'c 'Bx)
       0.0
       (symbolic-diff 'c 'Bx)
       0.0
      '(+ (* 0.0 c) (* c 0.0))
      (symbolic-diff 'Bx 'Bx)
      1.0
     '(+ (* (+ (* 0.0 c) (* c 0.0)) Bx) (* (* c c) 1.0))
    '(+ (* 0.0 (* (* c c) Bx)) (* b_fact (+ (* (+ (* 0.0 c) (* c 0.0)) Bx) (* (* c c) 1.0))))
    (symbolic-simp '(+ (* 0.0 (* (* c c) Bx)) (* b_fact (+ (* (+ (* 0.0 c) (* c 0.0)) Bx) (* (* c c) 1.0)))))
     (symbolic-simp-rule '(+ (* 0.0 (* (* c c) Bx)) (* b_fact (+ (* (+ (* 0.0 c) (* c 0.0)) Bx) (* (* c c) 1.0)))))
      (symbolic-simp-rule '(* 0.0 (* (* c c) Bx)))
      0.0
      (symbolic-simp-rule '(* b_fact (+ (* (+ (* 0.0 c) (* c 0.0)) Bx) (* (* c c) 1.0))))
      '(+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 1.0)))
     '(+ 0.0 (+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 1.0))))
    (symbolic-simp '(+ 0.0 (+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 1.0)))))
     (symbolic-simp-rule '(+ 0.0 (+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 1.0)))))
     '(+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 1.0)))
    (symbolic-simp '(+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 1.0))))
     (symbolic-simp-rule '(+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 1.0))))
      (symbolic-simp-rule '(* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)))
       (symbolic-simp-rule 'b_fact)
       'b_fact
       (symbolic-simp-rule '(* (+ (* 0.0 c) (* c 0.0)) Bx))
        (symbolic-simp-rule '(+ (* 0.0 c) (* c 0.0)))
         (symbolic-simp-rule '(* 0.0 c))
         0.0
         (symbolic-simp-rule '(* c 0.0))
         '(* 0.0 c)
        '(+ 0.0 (* 0.0 c))
        (symbolic-simp-rule 'Bx)
        'Bx
       '(* (+ 0.0 (* 0.0 c)) Bx)
      '(* b_fact (* (+ 0.0 (* 0.0 c)) Bx))
      (symbolic-simp-rule '(* b_fact (* (* c c) 1.0)))
       (symbolic-simp-rule 'b_fact)
       'b_fact
       (symbolic-simp-rule '(* (* c c) 1.0))
       '(* c (* c 1.0))
      '(* b_fact (* c (* c 1.0)))
     '(+ (* b_fact (* (+ 0.0 (* 0.0 c)) Bx)) (* b_fact (* c (* c 1.0))))
    (symbolic-simp '(+ (* b_fact (* (+ 0.0 (* 0.0 c)) Bx)) (* b_fact (* c (* c 1.0)))))
     (symbolic-simp-rule '(+ (* b_fact (* (+ 0.0 (* 0.0 c)) Bx)) (* b_fact (* c (* c 1.0)))))
      (symbolic-simp-rule '(* b_fact (* (+ 0.0 (* 0.0 c)) Bx)))
       (symbolic-simp-rule 'b_fact)
       'b_fact
       (symbolic-simp-rule '(* (+ 0.0 (* 0.0 c)) Bx))
        (symbolic-simp-rule '(+ 0.0 (* 0.0 c)))
        '(* 0.0 c)
        (symbolic-simp-rule 'Bx)
        'Bx
       '(* (* 0.0 c) Bx)
      '(* b_fact (* (* 0.0 c) Bx))
      (symbolic-simp-rule '(* b_fact (* c (* c 1.0))))
       (symbolic-simp-rule 'b_fact)
       'b_fact
       (symbolic-simp-rule '(* c (* c 1.0)))
        (symbolic-simp-rule 'c)
        'c
        (symbolic-simp-rule '(* c 1.0))
        '(* 1.0 c)
       '(* c (* 1.0 c))
      '(* b_fact (* c (* 1.0 c)))
     '(+ (* b_fact (* (* 0.0 c) Bx)) (* b_fact (* c (* 1.0 c))))
    (symbolic-simp '(+ (* b_fact (* (* 0.0 c) Bx)) (* b_fact (* c (* 1.0 c)))))
     (symbolic-simp-rule '(+ (* b_fact (* (* 0.0 c) Bx)) (* b_fact (* c (* 1.0 c)))))
      (symbolic-simp-rule '(* b_fact (* (* 0.0 c) Bx)))
       (symbolic-simp-rule 'b_fact)
       'b_fact
       (symbolic-simp-rule '(* (* 0.0 c) Bx))
       '(* 0.0 (* c Bx))
      '(* b_fact (* 0.0 (* c Bx)))
      (symbolic-simp-rule '(* b_fact (* c (* 1.0 c))))
       (symbolic-simp-rule 'b_fact)
       'b_fact
       (symbolic-simp-rule '(* c (* 1.0 c)))
       '(* 1.0 (* c c))
      '(* b_fact (* 1.0 (* c c)))
     '(+ (* b_fact (* 0.0 (* c Bx))) (* b_fact (* 1.0 (* c c))))
    (symbolic-simp '(+ (* b_fact (* 0.0 (* c Bx))) (* b_fact (* 1.0 (* c c)))))
     (symbolic-simp-rule '(+ (* b_fact (* 0.0 (* c Bx))) (* b_fact (* 1.0 (* c c)))))
      (symbolic-simp-rule '(* b_fact (* 0.0 (* c Bx))))
      '(* 0.0 (* b_fact (* c Bx)))
      (symbolic-simp-rule '(* b_fact (* 1.0 (* c c))))
      '(* 1.0 (* b_fact (* c c)))
     '(+ (* 0.0 (* b_fact (* c Bx))) (* 1.0 (* b_fact (* c c))))
    (symbolic-simp '(+ (* 0.0 (* b_fact (* c Bx))) (* 1.0 (* b_fact (* c c)))))
     (symbolic-simp-rule '(+ (* 0.0 (* b_fact (* c Bx))) (* 1.0 (* b_fact (* c c)))))
      (symbolic-simp-rule '(* 0.0 (* b_fact (* c Bx))))
      0.0
      (symbolic-simp-rule '(* 1.0 (* b_fact (* c c))))
      '(* b_fact (* c c))
     '(+ 0.0 (* b_fact (* c c)))
    (symbolic-simp '(+ 0.0 (* b_fact (* c c))))
     (symbolic-simp-rule '(+ 0.0 (* b_fact (* c c))))
     '(* b_fact (* c c))
    (symbolic-simp '(* b_fact (* c c)))
     (symbolic-simp-rule '(* b_fact (* c c)))
      (symbolic-simp-rule 'b_fact)
      'b_fact
      (symbolic-simp-rule '(* c c))
       (symbolic-simp-rule 'c)
       'c
       (symbolic-simp-rule 'c)
       'c
      '(* c c)
     '(* b_fact (* c c))
    '(* b_fact (* c c))
    (symbolic-diff '(* b_fact (* (* c c) Bx)) 'psi)
     (symbolic-diff 'b_fact 'psi)
     0.0
     (symbolic-diff '(* (* c c) Bx) 'psi)
      (symbolic-diff '(* c c) 'psi)
       (symbolic-diff 'c 'psi)
       0.0
       (symbolic-diff 'c 'psi)
       0.0
      '(+ (* 0.0 c) (* c 0.0))
      (symbolic-diff 'Bx 'psi)
      0.0
     '(+ (* (+ (* 0.0 c) (* c 0.0)) Bx) (* (* c c) 0.0))
    '(+ (* 0.0 (* (* c c) Bx)) (* b_fact (+ (* (+ (* 0.0 c) (* c 0.0)) Bx) (* (* c c) 0.0))))
    (symbolic-simp '(+ (* 0.0 (* (* c c) Bx)) (* b_fact (+ (* (+ (* 0.0 c) (* c 0.0)) Bx) (* (* c c) 0.0)))))
     (symbolic-simp-rule '(+ (* 0.0 (* (* c c) Bx)) (* b_fact (+ (* (+ (* 0.0 c) (* c 0.0)) Bx) (* (* c c) 0.0)))))
      (symbolic-simp-rule '(* 0.0 (* (* c c) Bx)))
      0.0
      (symbolic-simp-rule '(* b_fact (+ (* (+ (* 0.0 c) (* c 0.0)) Bx) (* (* c c) 0.0))))
      '(+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 0.0)))
     '(+ 0.0 (+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 0.0))))
    (symbolic-simp '(+ 0.0 (+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 0.0)))))
     (symbolic-simp-rule '(+ 0.0 (+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 0.0)))))
     '(+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 0.0)))
    (symbolic-simp '(+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 0.0))))
     (symbolic-simp-rule '(+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 0.0))))
      (symbolic-simp-rule '(* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)))
       (symbolic-simp-rule 'b_fact)
       'b_fact
       (symbolic-simp-rule '(* (+ (* 0.0 c) (* c 0.0)) Bx))
        (symbolic-simp-rule '(+ (* 0.0 c) (* c 0.0)))
         (symbolic-simp-rule '(* 0.0 c))
         0.0
         (symbolic-simp-rule '(* c 0.0))
         '(* 0.0 c)
        '(+ 0.0 (* 0.0 c))
        (symbolic-simp-rule 'Bx)
        'Bx
       '(* (+ 0.0 (* 0.0 c)) Bx)
      '(* b_fact (* (+ 0.0 (* 0.0 c)) Bx))
      (symbolic-simp-rule '(* b_fact (* (* c c) 0.0)))
       (symbolic-simp-rule 'b_fact)
       'b_fact
       (symbolic-simp-rule '(* (* c c) 0.0))
       '(* c (* c 0.0))
      '(* b_fact (* c (* c 0.0)))
     '(+ (* b_fact (* (+ 0.0 (* 0.0 c)) Bx)) (* b_fact (* c (* c 0.0))))
    (symbolic-simp '(+ (* b_fact (* (+ 0.0 (* 0.0 c)) Bx)) (* b_fact (* c (* c 0.0)))))
     (symbolic-simp-rule '(+ (* b_fact (* (+ 0.0 (* 0.0 c)) Bx)) (* b_fact (* c (* c 0.0)))))
      (symbolic-simp-rule '(* b_fact (* (+ 0.0 (* 0.0 c)) Bx)))
       (symbolic-simp-rule 'b_fact)
       'b_fact
       (symbolic-simp-rule '(* (+ 0.0 (* 0.0 c)) Bx))
        (symbolic-simp-rule '(+ 0.0 (* 0.0 c)))
        '(* 0.0 c)
        (symbolic-simp-rule 'Bx)
        'Bx
       '(* (* 0.0 c) Bx)
      '(* b_fact (* (* 0.0 c) Bx))
      (symbolic-simp-rule '(* b_fact (* c (* c 0.0))))
       (symbolic-simp-rule 'b_fact)
       'b_fact
       (symbolic-simp-rule '(* c (* c 0.0)))
        (symbolic-simp-rule 'c)
        'c
        (symbolic-simp-rule '(* c 0.0))
        '(* 0.0 c)
       '(* c (* 0.0 c))
      '(* b_fact (* c (* 0.0 c)))
     '(+ (* b_fact (* (* 0.0 c) Bx)) (* b_fact (* c (* 0.0 c))))
    (symbolic-simp '(+ (* b_fact (* (* 0.0 c) Bx)) (* b_fact (* c (* 0.0 c)))))
     (symbolic-simp-rule '(+ (* b_fact (* (* 0.0 c) Bx)) (* b_fact (* c (* 0.0 c)))))
      (symbolic-simp-rule '(* b_fact (* (* 0.0 c) Bx)))
       (symbolic-simp-rule 'b_fact)
       'b_fact
       (symbolic-simp-rule '(* (* 0.0 c) Bx))
       '(* 0.0 (* c Bx))
      '(* b_fact (* 0.0 (* c Bx)))
      (symbolic-simp-rule '(* b_fact (* c (* 0.0 c))))
       (symbolic-simp-rule 'b_fact)
       'b_fact
       (symbolic-simp-rule '(* c (* 0.0 c)))
       '(* 0.0 (* c c))
      '(* b_fact (* 0.0 (* c c)))
     '(+ (* b_fact (* 0.0 (* c Bx))) (* b_fact (* 0.0 (* c c))))
    (symbolic-simp '(+ (* b_fact (* 0.0 (* c Bx))) (* b_fact (* 0.0 (* c c)))))
     (symbolic-simp-rule '(+ (* b_fact (* 0.0 (* c Bx))) (* b_fact (* 0.0 (* c c)))))
      (symbolic-simp-rule '(* b_fact (* 0.0 (* c Bx))))
      '(* 0.0 (* b_fact (* c Bx)))
      (symbolic-simp-rule '(* b_fact (* 0.0 (* c c))))
      '(* 0.0 (* b_fact (* c c)))
     '(+ (* 0.0 (* b_fact (* c Bx))) (* 0.0 (* b_fact (* c c))))
    (symbolic-simp '(+ (* 0.0 (* b_fact (* c Bx))) (* 0.0 (* b_fact (* c c)))))
     (symbolic-simp-rule '(+ (* 0.0 (* b_fact (* c Bx))) (* 0.0 (* b_fact (* c c)))))
      (symbolic-simp-rule '(* 0.0 (* b_fact (* c Bx))))
      0.0
      (symbolic-simp-rule '(* 0.0 (* b_fact (* c c))))
      0.0
     '(+ 0.0 0.0)
    (symbolic-simp '(+ 0.0 0.0))
     (symbolic-simp-rule '(+ 0.0 0.0))
     0.0
    (symbolic-simp 0.0)
     (symbolic-simp-rule 0.0)
     0.0
    0.0
   '((* b_fact (* c c)) 0.0)
  (symbolic-jacobian '((* b_fact (* c c)) 0.0) '(Bx psi))
   (symbolic-diff '(* b_fact (* c c)) 'Bx)
    (symbolic-diff 'b_fact 'Bx)
    0.0
    (symbolic-diff '(* c c) 'Bx)
     (symbolic-diff 'c 'Bx)
     0.0
     (symbolic-diff 'c 'Bx)
     0.0
    '(+ (* 0.0 c) (* c 0.0))
   '(+ (* 0.0 (* c c)) (* b_fact (+ (* 0.0 c) (* c 0.0))))
   (symbolic-simp '(+ (* 0.0 (* c c)) (* b_fact (+ (* 0.0 c) (* c 0.0)))))
    (symbolic-simp-rule '(+ (* 0.0 (* c c)) (* b_fact (+ (* 0.0 c) (* c 0.0)))))
     (symbolic-simp-rule '(* 0.0 (* c c)))
     0.0
     (symbolic-simp-rule '(* b_fact (+ (* 0.0 c) (* c 0.0))))
     '(+ (* b_fact (* 0.0 c)) (* b_fact (* c 0.0)))
    '(+ 0.0 (+ (* b_fact (* 0.0 c)) (* b_fact (* c 0.0))))
   (symbolic-simp '(+ 0.0 (+ (* b_fact (* 0.0 c)) (* b_fact (* c 0.0)))))
    (symbolic-simp-rule '(+ 0.0 (+ (* b_fact (* 0.0 c)) (* b_fact (* c 0.0)))))
    '(+ (* b_fact (* 0.0 c)) (* b_fact (* c 0.0)))
   (symbolic-simp '(+ (* b_fact (* 0.0 c)) (* b_fact (* c 0.0))))
    (symbolic-simp-rule '(+ (* b_fact (* 0.0 c)) (* b_fact (* c 0.0))))
     (symbolic-simp-rule '(* b_fact (* 0.0 c)))
     '(* 0.0 (* b_fact c))
     (symbolic-simp-rule '(* b_fact (* c 0.0)))
      (symbolic-simp-rule 'b_fact)
      'b_fact
      (symbolic-simp-rule '(* c 0.0))
      '(* 0.0 c)
     '(* b_fact (* 0.0 c))
    '(+ (* 0.0 (* b_fact c)) (* b_fact (* 0.0 c)))
   (symbolic-simp '(+ (* 0.0 (* b_fact c)) (* b_fact (* 0.0 c))))
    (symbolic-simp-rule '(+ (* 0.0 (* b_fact c)) (* b_fact (* 0.0 c))))
     (symbolic-simp-rule '(* 0.0 (* b_fact c)))
     0.0
     (symbolic-simp-rule '(* b_fact (* 0.0 c)))
     '(* 0.0 (* b_fact c))
    '(+ 0.0 (* 0.0 (* b_fact c)))
   (symbolic-simp '(+ 0.0 (* 0.0 (* b_fact c))))
    (symbolic-simp-rule '(+ 0.0 (* 0.0 (* b_fact c))))
    '(* 0.0 (* b_fact c))
   (symbolic-simp '(* 0.0 (* b_fact c)))
    (symbolic-simp-rule '(* 0.0 (* b_fact c)))
    0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff '(* b_fact (* c c)) 'psi)
    (symbolic-diff 'b_fact 'psi)
    0.0
    (symbolic-diff '(* c c) 'psi)
     (symbolic-diff 'c 'psi)
     0.0
     (symbolic-diff 'c 'psi)
     0.0
    '(+ (* 0.0 c) (* c 0.0))
   '(+ (* 0.0 (* c c)) (* b_fact (+ (* 0.0 c) (* c 0.0))))
   (symbolic-simp '(+ (* 0.0 (* c c)) (* b_fact (+ (* 0.0 c) (* c 0.0)))))
    (symbolic-simp-rule '(+ (* 0.0 (* c c)) (* b_fact (+ (* 0.0 c) (* c 0.0)))))
     (symbolic-simp-rule '(* 0.0 (* c c)))
     0.0
     (symbolic-simp-rule '(* b_fact (+ (* 0.0 c) (* c 0.0))))
     '(+ (* b_fact (* 0.0 c)) (* b_fact (* c 0.0)))
    '(+ 0.0 (+ (* b_fact (* 0.0 c)) (* b_fact (* c 0.0))))
   (symbolic-simp '(+ 0.0 (+ (* b_fact (* 0.0 c)) (* b_fact (* c 0.0)))))
    (symbolic-simp-rule '(+ 0.0 (+ (* b_fact (* 0.0 c)) (* b_fact (* c 0.0)))))
    '(+ (* b_fact (* 0.0 c)) (* b_fact (* c 0.0)))
   (symbolic-simp '(+ (* b_fact (* 0.0 c)) (* b_fact (* c 0.0))))
    (symbolic-simp-rule '(+ (* b_fact (* 0.0 c)) (* b_fact (* c 0.0))))
     (symbolic-simp-rule '(* b_fact (* 0.0 c)))
     '(* 0.0 (* b_fact c))
     (symbolic-simp-rule '(* b_fact (* c 0.0)))
      (symbolic-simp-rule 'b_fact)
      'b_fact
      (symbolic-simp-rule '(* c 0.0))
      '(* 0.0 c)
     '(* b_fact (* 0.0 c))
    '(+ (* 0.0 (* b_fact c)) (* b_fact (* 0.0 c)))
   (symbolic-simp '(+ (* 0.0 (* b_fact c)) (* b_fact (* 0.0 c))))
    (symbolic-simp-rule '(+ (* 0.0 (* b_fact c)) (* b_fact (* 0.0 c))))
     (symbolic-simp-rule '(* 0.0 (* b_fact c)))
     0.0
     (symbolic-simp-rule '(* b_fact (* 0.0 c)))
     '(* 0.0 (* b_fact c))
    '(+ 0.0 (* 0.0 (* b_fact c)))
   (symbolic-simp '(+ 0.0 (* 0.0 (* b_fact c))))
    (symbolic-simp-rule '(+ 0.0 (* 0.0 (* b_fact c))))
    '(* 0.0 (* b_fact c))
   (symbolic-simp '(* 0.0 (* b_fact c)))
    (symbolic-simp-rule '(* 0.0 (* b_fact c)))
    0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff 0.0 'Bx)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff 0.0 'psi)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
  '((0.0 0.0) (0.0 0.0))
  (symbolic-eigvals2 '((0.0 0.0) (0.0 0.0)))
  '(0.0 0.0)
  (symbolic-eigvals2 '((0.0 0.0) (0.0 0.0)))
  '(0.0 0.0)
  (symbolic-simp 0.0)
   (symbolic-simp-rule 0.0)
   0.0
  0.0
  (symbolic-simp 0.0)
   (symbolic-simp-rule 0.0)
   0.0
  0.0
  (symbolic-simp 0.0)
   (symbolic-simp-rule 0.0)
   0.0
  0.0
  (symbolic-simp 0.0)
   (symbolic-simp-rule 0.0)
   0.0
  0.0
  (is-real 1.0 '(Bx psi) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-real 1.0 '(Bx psi) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-real 1.0 '(Bx psi) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-real 0.0 '(Bx psi) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-real '(cond ((< x 0.0) 0.5) (else -0.5)) '(Bx psi) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
   (is-real 0.5 '(Bx psi) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
   #t
  (is-real -0.5 '(Bx psi) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-non-negative 0.0 '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-non-negative 0.0 '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-non-negative 0.0 '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-non-negative 0.0 '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
 #t
