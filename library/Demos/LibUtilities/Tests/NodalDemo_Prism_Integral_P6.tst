<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Nodal prism integration, evenly spaced points, P = 6</description>
    <executable>NodalDemo</executable>
    <parameters>--order 6 --type 27 --integral</parameters>
    <metrics>
        <metric type="Linf" id="1">
            <value tolerance="1e-12">3.00544e-06</value>
        </metric>
    </metrics>
</test>
