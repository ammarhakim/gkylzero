// Exercise cluster startup against a real reporting repository without the new helper.
// Jenkins steps are stubbed; Git extraction and loaded helper configuration are real.
def root = new File('ci/jenkins').canonicalFile
def temp = File.createTempDir('cluster-reporting-pin-', '')
def run = { String command, File cwd, Map environment ->
    def process = new ProcessBuilder('bash', '-c', command).directory(cwd)
    environment.each { key, value -> process.environment().put(key.toString(), value.toString()) }
    def child = process.start()
    def out = new StringBuffer(), err = new StringBuffer()
    child.waitForProcessOutput(out, err)
    assert child.exitValue() == 0: err.toString()
    out.toString()
}
try {
    def repo = new File(temp, 'reporting-repo')
    def files = new File(repo, 'ci/jenkins')
    files.mkdirs()
    new File(files, 'jenkins_reporting.groovy').text = '''
        def configure(Map settings) { binding.setVariable('settings', settings) }
        return this
    '''
    ['github_report.py', 'baseline_cache.sh', 'prebuilt_config.py'].each {
        new File(files, it).text = '# old trusted tool\n'
    }
    run('git init -q && git add . && git -c user.name=Fixture -c user.email=fixture@example.invalid commit -qm fixture', repo, [:])
    def pin = run('git rev-parse HEAD', repo, [:]).trim()
    ['stellar_cpu', 'perlmutter_gpu'].each { platform ->
        def source = new File(root, "jenkinsfile.${platform}").text
        source = source.substring(0, source.indexOf('// Manual, CPU-only') >= 0
            ? source.indexOf('// Manual, CPU-only') : source.indexOf('def buildsToKeep'))
        def workspace = new File(temp, platform)
        workspace.mkdirs()
        def cwd = workspace
        def environment = [GKEYLL_CI_TRUSTED_REF: pin]
        def trustedScm = [trusted: true]
        def checkouts = [], loads = []
        def binding = new Binding([
            env: environment, scm: trustedScm,
            pwd: { Map options = [:] -> options.tmp ? workspace.path + '@tmp' : cwd.path },
            dir: { String path, Closure body ->
                def previous = cwd
                cwd = path.startsWith('/') ? new File(path) : new File(cwd, path)
                cwd.mkdirs()
                try { body() } finally { cwd = previous }
            },
            deleteDir: { cwd.deleteDir(); cwd.mkdirs() },
            withEnv: { List entries, Closure body ->
                def previous = new HashMap(environment)
                entries.each { def pair = it.toString().split('=', 2); environment[pair[0]] = pair[1] }
                try { body() } finally { environment.clear(); environment.putAll(previous) }
            },
            sh: { Object options ->
                def command = options instanceof Map ? options.script : options
                run(command.toString().replace('https://github.com/gkeyllorg/gkeyll.git', repo.path), cwd, environment)
            },
            checkout: { Object spec ->
                assert spec.is(trustedScm)
                assert cwd.path.startsWith(workspace.path + '@tmp/')
                checkouts << spec
                def target = new File(cwd, 'ci/jenkins/jenkins_pipeline.groovy')
                target.parentFile.mkdirs()
                target.text = new File(root, 'jenkins_pipeline.groovy').text
            },
            load: { String path ->
                loads << path
                new GroovyShell().evaluate(new File(cwd, path))
            },
            error: { String message -> throw new IllegalArgumentException(message) }
        ])
        def pipeline = new GroovyShell(binding).parse(source)
        def settings = [:]
        pipeline.loadReporting(settings)
        assert settings.reportingCommit == pin
        assert checkouts.size() == 1
        assert loads == ['jenkins_reporting.groovy', 'ci/jenkins/jenkins_pipeline.groovy']
        assert pipeline.pipelineHelpers.reporting.is(pipeline.reporting)
    }
} finally { temp.deleteDir() }
println('Cluster startup succeeds with a reporting pin predating shared Pipeline helpers')
