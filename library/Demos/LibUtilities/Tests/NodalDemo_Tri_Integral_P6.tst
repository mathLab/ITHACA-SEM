<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Nodal triangle integration, Fekete points, P = 6</description>
    <executable>NodalDemo</executable>
    <parameters>--order 6 --type 23 --integral</parameters>
    <metrics>
        <metric type="Linf" id="1">
            <value tolerance="1e-12">7.28803e-07</value>
        </metric>
    </metrics>
</test>
