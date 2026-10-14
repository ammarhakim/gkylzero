// Offline contracts for the shared helper and its machine entrypoints.
def root = new File('ci/jenkins')
def calls = []
def environment = [WORKSPACE: '/fixture']
def binding = new Binding([
    env: environment,
    checkout: { Map spec -> calls << spec },
    sh: { Object command -> calls << command; 'a' * 40 + '\n' },
    error: { String message -> throw new IllegalArgumentException(message) },
    withEnv: { List values, Closure body -> calls << values; body() }
])
def helper = new GroovyShell(binding).parse(new File(root, 'jenkins_pipeline.groovy'))
helper.configure([loggedSh: { String label, String command -> calls << [label, command] }])
helper.checkoutRef('feature/test', 'candidate')
assert calls[-1].userRemoteConfigs[0].refspec == '+refs/heads/feature/test:refs/remotes/origin/candidate'
helper.checkoutRef('b' * 40, 'baseline')
assert calls[-1].userRemoteConfigs[0].refspec == '+' + 'b' * 40 + ':refs/remotes/origin/baseline'
helper.checkoutPullRequest('42')
assert calls[-1].userRemoteConfigs[0].refspec == '+refs/pull/42/head:refs/remotes/origin/candidate'
assert helper.resolveRefCommit('A' * 40) == 'a' * 40
assert helper.resolveRefCommit('main') == 'a' * 40
helper.configureAndCompileCRegressions('candidate', '/gkeyll', '/source', '/prefix', '/modules', '8', '--parallel', '--mpiexec mpiexec')
assert calls[-2].collect { it.toString() }.contains('CI_REGRESSION_COMPILE_JOBS=8')
assert calls[-2].collect { it.toString() }.contains('CI_REGRESSION_RUN_OPTIONS=--parallel')
assert calls[-2].collect { it.toString() }.contains('CI_REGRESSION_CONFIGURE_OPTIONS=--mpiexec mpiexec')
assert calls[-1][0] == 'candidate-c-compile'
helper.finalizeUnfinishedSlurmJob()
assert calls[-1].contains('scancel "$job_id"')
['personal', 'stellar_cpu', 'perlmutter_gpu'].each { platform ->
    def source = new File(root, "jenkinsfile.${platform}").text
    new GroovyShell().parse(source)
    assert source.contains('git show FETCH_HEAD:ci/jenkins/jenkins_pipeline.groovy > jenkins_pipeline.groovy')
    assert source.contains("pipelineHelpers = load('jenkins_pipeline.groovy')")
    assert source.contains('pipelineHelpers.configure(reporting)')
    assert source.contains('return pipelineHelpers.usePrebuiltConfig(prefix)')
}
println('Shared Pipeline helper tests passed')
