FROM sphinxdoc/sphinx as builder

WORKDIR /docs
COPY ./ ./
RUN rm -rf docs/build

RUN pip install sphinx_rtd_theme
RUN pip install .

RUN cd docs && make html

FROM nginx:alpine
COPY --from=builder /docs/docs/build/html /usr/share/nginx/html