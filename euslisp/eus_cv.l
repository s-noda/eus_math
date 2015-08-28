#!/usr/bin/env roseus

(defvar *cv-plugin*
  (labels
   ((library_search
     (str &key
	  (depth 0)
	  colon-pos lib-path)
     (format t "  [~A] target=" depth)
     (cond
      ((eq (length str) 0)
       (format t "~% euscv_plugin_not_found in eus_cv.l~%")
       (exit -1))
      ((and (setq colon-pos (or (position #\: str) (length str)))
	    (setq lib-path (subseq str 0 colon-pos))
	    (setq lib-path
		  (if (eq (aref lib-path (- (length lib-path) 1)) #\/)
		      (subseq lib-path 0 (- (length lib-path) 1))
		    lib-path))
	    (probe-file (setq lib-path
			      (print (format nil "~A/libeus_cv.so" lib-path)))))
       (load-foreign lib-path))
      (t
       (library_search (subseq str (min (length str) (+ colon-pos 1)))
		       :depth (+ depth 1))))))
   (library_search (format nil "~A:~A/lib"
			   (unix:getenv "LD_LIBRARY_PATH")
			   (read-line (piped-fork "rospack find eus_math"))))))

(defforeign _cv-kmeans
  *cv-plugin*
  "cv_kmeans"
  (:integer :integer :integer :string :string :string :integer :double :integer)
  :integer)

(defun cv-kmeans
  (&key
   (k 3)
   (sample-n 10)
   (length 3)
   (points (instantiate float-vector (* sample-n length)))
   (id (instantiate integer-vector sample-n))
   (centers (instantiate float-vector (* k length)))
   (max-iter 10)
   (eps 0.01)
   ;;
   (points-str (dvector2float-bytestring points))
   (id-str (lvector2integer-bytestring id))
   (centers-str (dvector2float-bytestring centers))
   ;;
   (debug? t)
   )
  (_cv-kmeans k sample-n length points-str id-str centers-str max-iter eps (if debug? 1 0))
  ;;
  (if points (float-bytestring2dvector points-str points))
  (if centers (float-bytestring2dvector centers-str centers))
  (if id (integer-bytestring2lvector id-str id))
  ;;
  ;; (print points) (print centers) (print id)
  )

(defun test-cv-kmeans
  (&key
   (k 3)
   (sample-n 10)
   (length 3)
   (points (let* ((p (instantiate float-vector (* sample-n length)))
		  (center (instantiate float-vector length)))
	     (dotimes (j length)
	       (dotimes (i sample-n)
		 (setf (aref p (+ i (* j sample-n)))
		       (+ (aref center j) (random 1.0)))
		 )
	       (dotimes (j (length center))
		 (setf (aref center j) (+ (aref center j) 100))))
	     p))
   (id (instantiate integer-vector sample-n))
   (centers (instantiate float-vector (* k length)))
   (max-iter 10)
   (eps 0.01)
   )
  (cv-kmeans :k k :sample-n sample-n :length length :points points :id id :centers centers :max-iter max-iter :eps eps)
  (print points)
  (print id)
  (print centers))

