(ns propeller.gp
  (:require [clojure.string]
            [clojure.pprint]
            [propeller.genome :as genome]
            [propeller.variation :as variation]
            [propeller.push.instructions.bool]
            [propeller.push.instructions.character]
            [propeller.push.instructions.code]
            [propeller.push.instructions.input-output]
            [propeller.push.instructions.numeric]
            [propeller.push.instructions.polymorphic]
            [propeller.push.instructions.string]
            [propeller.push.instructions.vector]))

(defn report
  "Reports information each generation."
  [pop generation argmap]
  (let [best (first pop)]
    (clojure.pprint/pprint {:generation            generation
                            :best-plushy           (:plushy best)
                            :best-program          (genome/plushy->push (:plushy best) argmap)
                            :best-total-error      (:total-error best)
                            :best-errors           (:errors best)
                            :best-behaviors        (:behaviors best)
                            :genotypic-diversity   (float (/ (count (distinct (map :plushy pop))) (count pop)))
                            :behavioral-diversity  (float (/ (count (distinct (map :behaviors pop))) (count pop)))
                            :average-genome-length (float (/ (reduce + (map count (map :plushy pop))) (count pop)))
                            :average-total-error   (float (/ (reduce + (map :total-error pop)) (count pop)))})
    (println)))

(defn track-solutions
  "Track all solutions we've found so far.
   struct-solutions is a map where the keys are the hash values of solution
   programs, and the values are how many times we've seen them so far."
  [inputs population struct-solutions]
  (let [solutions (filter #(zero? (:total-error %)) population)]
    (loop [solutions solutions
           structs struct-solutions]
      (if (empty? solutions)
        structs
        (let [push-prog (genome/plushy->push (:plushy (first solutions)))
              prog-hash (hash push-prog)
              behavior ((partial error-function argmap 
                                 inputs) (first solutions))
              
              ]
          (recur (rest solutions)
                 (if (get structs prog-hash)
                   (update structs prog-hash inc)
                   (assoc structs prog-hash 1))))))))

(def seed-solution
  '(:in1 :in2 :in3 :in4 :integer_min :integer_min :integer_min :integer_print))

(defn gp
  "Main GP loop."
  [{:keys [population-size max-generations error-function instructions
           max-initial-plushy-size random-inputs-for-generalizability]
    :as   argmap}]
  ;;
  (prn {:starting-args (update (update argmap :error-function str) :instructions str)})
  (println)
  ;;
  (loop [generation 0
         population (conj
                     (repeatedly
                           (dec population-size)
                           #(hash-map :plushy (genome/make-random-plushy
                                               instructions
                                               max-initial-plushy-size)))
                     {:plushy seed-solution})
         struct-solutions {}]
    (let [evaluated-pop (sort-by :total-error
                                 (#?(:clj  pmap
                                     :cljs map)
                                   (partial error-function argmap (:training-data argmap))
                                   population))
          best-individual (first evaluated-pop)
          struct-solutions (track-solutions random-inputs-for-generalizability
                                            evaluated-pop 
                                            struct-solutions)]
      (if (:custom-report argmap)
        ((:custom-report argmap) evaluated-pop generation argmap)
        (report evaluated-pop generation argmap))
      (prn "Struct Solutions Map:" struct-solutions)
      (cond
        ;; Success on training cases is verified on testing cases
         ;(zero? (:total-error best-individual))
        false ;; This is to make it so we don't stop when finding a solution
        (do (prn {:success-generation generation})
            (prn {:total-test-error
                  (:total-error (error-function argmap (:testing-data argmap) best-individual))})
            (#?(:clj shutdown-agents)))
        ;;
        (>= generation max-generations)
        nil
        ;;
        :else (recur (inc generation)
                     (if (:elitism argmap)
                       (conj (repeatedly (dec population-size)
                                         #(variation/new-individual evaluated-pop argmap))
                             (first evaluated-pop))
                       (repeatedly population-size
                                   #(variation/new-individual evaluated-pop argmap)))
                     struct-solutions)))))
