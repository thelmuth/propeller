(ns propeller.problems.PSB2.luhn
  (:require [psb2.core :as psb2]
            [propeller.genome :as genome]
            [propeller.push.interpreter :as interpreter]
            [propeller.utils :as utils]
            [propeller.push.instructions :refer [get-stack-instructions]]
            [propeller.push.state :as state]
            [propeller.tools.math :as math]
            [propeller.gp :as gp]
            #?(:cljs [cljs.reader :refer [read-string]])))

; ===========  PROBLEM DESCRIPTION  ============================
; LUHN from PSB2
; Given a vector of 16 digits, implement Luhn’s
; algorithm to verify a credit card number, such that it follows
; the following rules: double every other digit starting with
; the second digit. If any of the results are over 9, subtract 9
; from them. Return the sum of all of the new digits.
;
; Source: https://arxiv.org/pdf/2106.06086.pdf
; ===============================================================

(def train-and-test-data (psb2/fetch-examples "data" "luhn" 200 2000))

; Random integer between -100 and 100 (from smallest)
(defn random-int [] (- (rand-int 201) 100))

(def instructions
  (utils/not-lazy
    (concat
      ;;; stack-specific instructions
      (get-stack-instructions #{:exec :integer :boolean :vector_integer :print})
      ;;; input instructions
      (list :in1)
      ;;; close
      (list 'close)
      ;;; ERCs (constants)
      (list 0 2 9 10 random-int))))

(defn error-function
  [argmap data individual]
  (let [program (genome/plushy->push (:plushy individual) argmap)
        inputs (map (fn [i] (get i :input1)) data)
        correct-outputs (map (fn [i] (get i :output1)) data)
        outputs (map (fn [input]
                       (state/peek-stack
                         (interpreter/interpret-program
                           program
                           (assoc state/empty-state :input {:in1 input})
                           (:step-limit argmap))
                         :integer))
                     inputs)
        errors (map (fn [correct-output output]
                      (if (= output :no-stack-item)
                        1000000
                        (math/abs (- correct-output output))))
                    correct-outputs
                    outputs)]
    (assoc individual
      :behaviors outputs
      :errors errors
      :total-error #?(:clj  (apply +' errors)
                      :cljs (apply + errors)))))

(defn -main
  "Runs propel-gp, giving it a map of arguments."
  [& args]
  (gp/gp
    (merge
      {:instructions            instructions
       :error-function          error-function
       :training-data           (:train train-and-test-data)
       :testing-data            (:test train-and-test-data)
       :max-generations         300
       :population-size         1000
       :max-initial-plushy-size 250
       :step-limit              2000
       :parent-selection        :lexicase
       :tournament-size         5
       :umad-rate               0.1
       :variation               {:umad 1.0 :crossover 0.0}
       :elitism                 false}
      (apply hash-map (map #(if (string? %) (read-string %) %) args))))
  (#?(:clj shutdown-agents)))
