<template>
  <q-page class="row justify-center wrap">
    <q-ajax-bar ref="bar" position="bottom" color="primary" size="10px" skip-hijack />
    <q-card flat bordered class="col-xs-12 col-md-10 bg-grey-1">
      <q-card-section>
        <div class="text-h3 text-bold tex-italic">LU decomposition</div>
      </q-card-section>

      <q-card-section class="flex" ref="description">
       <vue-mathjax ref="answer" :formula="description" :options="props" class="col-12" />
      </q-card-section>
      <div class="q-pa-md">
        <q-card class="sub-card bg-grey-4">
          <q-card-section>
            <div class="text-h6">Configuration</div>
          </q-card-section>

          <q-card-section>
            <div class="q-pa-md">
              <div class="q-gutter-md" style>
                <q-select
                  outlined
                  v-model="type"
                  label="Pick type of problem to solve"
                  :options="types"
                  :dense="dense"
                  :options-dense="denseOpts"
                >
                  <template v-slot:append>
                    <q-icon name="close" @click.stop="model = ''" class="cursor-pointer" />
                  </template>
                </q-select>

                <q-select
                  outlined
                  v-model="method"
                  label="Pick a method"
                  :options="methods"
                  :dense="dense"
                  :options-dense="denseOpts"
                >
                  <template v-slot:append>
                    <q-icon name="close" @click.stop="model = ''" class="cursor-pointer" />
                  </template>
                </q-select>

                <q-input
                  bottom-slots
                  v-model="n"
                  type="number"
                  label="Pick dimension"
                  counter
                  :dense="dense"
                >
                  <template v-slot:append>
                    <q-icon name="close" @click="n = ''" class="cursor-pointer" />
                  </template>

                  <template v-slot:hint>matrix dimension i.e 2 for 2 by 2 matrix</template>
                </q-input>
              </div>
            </div>
          </q-card-section>

          <q-separator dark />

          <q-card-actions>
            <q-btn :loading="loading" color="red" @click="simulateProgress">
              Generate Question
              <template v-slot:loading>
                <q-spinner-gears />
              </template>
            </q-btn>
          </q-card-actions>
        </q-card>
      </div>
      <div class="q-pa-md" v-show="question">
        <q-card class="sub-card bg-grey-4">
          <q-card-section>
            <div class="text-h6">Question</div>
          </q-card-section>

          <q-card-section>
            <vue-mathjax ref="question" :formula="question" :options="props" class="col-12" />
          </q-card-section>

          <q-separator dark />

          <q-card-actions>
            <q-toggle
              v-model="showAnswer"
              checked-icon="check"
              color="green"
              unchecked-icon="clear"
              label="show answer"
            />
            <q-toggle
              v-model="showWorking"
              checked-icon="check"
              color="red"
              label="show steps"
              unchecked-icon="clear"
            />
          </q-card-actions>
        </q-card>
      </div>
      <div class="q-pa-md" v-show="showAnswer">
        <q-card class="sub-card bg-grey-4">
          <q-card-section>
            <div class="text-h6">Answer</div>
          </q-card-section>

          <q-card-section>
            <vue-mathjax ref="answer" :formula="answer" :options="props" class="col-12" />
          </q-card-section>
        </q-card>
      </div>
      <div class="q-pa-md" v-show="showWorking">
        <q-card class="sub-card bg-grey-4">
          <q-card-section>
            <div class="text-h6">working</div>
          </q-card-section>

          <q-card-section>
            <vue-mathjax ref="working" :formula="working" :options="props" class="col-12" />
          </q-card-section>
        </q-card>
      </div>
    </q-card>
  </q-page>
</template>
<script>
import { VueMathjax } from "vue-mathjax";
import { mapState } from "vuex";
export default {
  name: "index",
  components: {
    VueMathjax
  },
  data() {
    return {
      answer: undefined,
      step: 1,
      loading: false,
      type: "Decompose A = LU",
      method: "DoLittle",
      description: " This solver is meant to help students in learning how to solve the linear system of equations $Ax = b$ using $LU$ factorisation",
      types: ["Decompose A = LU", "Solve AX = B"],
      n: 2,
      question: undefined,
      methods: ["DoLittle", "Crout", "Naive Gauss"],
      dense: false,
      denseOpts: false,
      props: {
        type: Array,
        default: () => {}
      },
      working: undefined,
      showWorking: false,
      showAnswer: false
    };
  },
  computed: {
    ...mapState("app", ["dialog", "maximizedToggle", "message"]),
    logoUrl() {
      return require("@/assets/isolated-monochrome-black.svg");
      // The path could be '../assets/img.png', etc., which depends on where your vue file is
    },
    dialog: {
      get() {
        return this.$store.state.app.dialog;
      },
      set(val) {
        this.$store.commit("app/setdialog", val);
      }
    },
    message: {
      get() {
        return this.$store.state.app.message;
      },
      set(val) {
        this.$store.commit("app/setMessage", val);
      }
    },
    maximizedToggle: {
      get() {
        return this.$store.state.app.maximizedToggle;
      },
      set(val) {
        this.$store.commit("app/setMaximizedToggle", val);
      }
    }
  },
  methods: {
    setMessage(message) {
      this.$store.commit("app/setMessage", message);
    },
    toggleDialog() {
      this.$store.commit("app/setDialog", !this.$store.state.app.dialog);
    },
    simulateProgress() {
      const bar = this.$refs.bar;
      bar.start();
      this.generateQuestion();
    },
    generateQuestion() {
      let payload = {
        m: `${this.n}`,
        method: `${this.methods.indexOf(this.method) + 1}`,
        type: `${this.types.indexOf(this.type) + 1}`
      };
      this.$store.dispatch("post/getQuestionSolution", payload).then(
        response => {
          this.question = response.data["question"];
          // this.update_element("question");
          this.answer = response.data["answer"]
          this.working = response.data["working"];
          // this.update_element("working");
          if (this.$refs.bar) {
            this.$refs.bar.stop();
          }
        },
        error => {
          this.setMessage(error.response.data);
          this.toggleDialog();
          this.reset();
        }
      );
    },
    reset() {
      this.question = undefined;
      this.working = undefined;
    },
    update_element(element) {
      this.$nextTick(() => {
        window.MathJax.Hub.Queue([
          "Typeset",
          window.MathJax.Hub,
          this.$refs[`${element}`]
        ]);
      });
    }
  }
};
</script>
<style lang="css" scoped>
.sub-card {
  overflow: auto;
}
</style>