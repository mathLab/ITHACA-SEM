<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Nodal prism interpolation, evenly spaced points, P = 7</description>
    <executable>NodalDemo</executable>
    <parameters>--order 7 --type 27 --interp -0.3,-0.636,-0.9</parameters>
    <metrics>
        <metric type="Linf" id="1">
            <value tolerance="1e-12">3.18397e-08</value>
        </metric>
    </metrics>
</test>
